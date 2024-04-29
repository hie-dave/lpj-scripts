from typing import Callable
from multiprocessing import BoundedSemaphore
from ozflux_logging import *
import concurrent.futures, mpi4py.MPI, mpi4py.futures
from ozflux_common import Task, Job
from traceback import format_exception

# Rank of the master node.
_RANK_MASTER = 0

# Tag used for transmission of commands from master to worker nodes.
_TAG_CMD = 0

# MPI tag used for communication of worker rank.
_TAG_RANK = 2

# MPI tag used for communication of progress reports.
_TAG_PROGRESS = 3

# MPI tag used by workers to signal that they are in an idle state.
_TAG_WORKER_IDLE = 4

# Tag used for transmission of job function from master to worker nodes.
_TAG_TASK = 5

# Tag used for transmission of job error details from worker to master.
_TAG_ERROR = 6

# The command send by the master node to signal that no more work needs to be
# done.
_CMD_DONE = 0

# Command received from the master node which indicates that a job is to be run.
_CMD_JOB = 1

# The communication world.
_comm = mpi4py.MPI.COMM_WORLD

def _mpi_worker_report_progress(progress: float):
	"""
	Write a progress report.

	@param progress: Current progress [0, 1].
	"""
	# Non-blocking send.
	_comm.isend(progress, _RANK_MASTER, tag = _TAG_PROGRESS)

def _mpi_master_send_job(worker: int, task: Task):
	"""
	Send a job to the specified worker.

	@param worker: Rank of the worker which will execute the job.
	@param task: The task to be executed.
	"""
	# Send the function to be executed.
	log_diagnostic(f"Sending CMD_JOB to worker {worker}")
	_comm.send(_CMD_JOB, worker, _TAG_CMD)
	log_diagnostic(f"Sending job to worker {worker}")
	_comm.send(task, worker, _TAG_TASK)
	log_diagnostic(f"Job successfully submitted to worker {worker}")

class MPIWorker:
	"""
	Stores metadata and handles communication with an MPI worker.
	"""
	def __init__(self, rank: int):
		self.rank = rank
		# Initialise _idle to false to ensure that we wait for an idle
		# notification to be received from the worker before we send it any
		# jobs. If we don't wait for such a request, we would deadlock.
		self._idle = False
		self.progress: float = 0.0
		self._idle_request: mpi4py.MPI.Request = None
		self._progress_request: mpi4py.MPI.Request = None
		self._error_request: mpi4py.MPI.Request = None
		self.job: Job = None

	def is_idle(self) -> bool:
		"""
		Check if this worker is idle.
		"""
		# If the worker is already idle, return true immediately.
		if self._idle:
			return True

		# Worker is currently busy. Create a non-blocking request for a
		# WORKER_IDLE message.
		if self._idle_request is None:
			log_diagnostic(f"Creating idle request for worker {self.rank}")
			self._idle_request = _comm.irecv(source = self.rank, tag = _TAG_WORKER_IDLE)
  
		# Ensure we also check for error notifications from the worker.
		if self._error_request is None:
			log_diagnostic(f"Creating error request for worker {self.rank}")
			self._error_request = _comm.irecv(source = self.rank, tag = _TAG_ERROR)

		(has_error, error) = self._error_request.test()
		if (has_error):
			if error is None:
				log_error("<<< PROTOCOL ERROR: received None exception from worker >>>")
			else:
				log_error(f"Worker {self.rank} has run with error:")
				log_error(error)
			self._error_request = None
			_comm.Abort(1)

		# If the request for a WORKER_IDLE message has completed, the worker is
		# now idle.
		(is_idle, _) = self._idle_request.test()
		if is_idle:
			log_diagnostic(f"Worker {self.rank} is idle")
			self._idle_request = None
			self._idle = True

		return self._idle

	def has_progress_update(self) -> bool:
		"""
		Check if a progress update has been received from this worker.
		"""
		if self._idle:
			return False
		if self._progress_request is None:
			self._progress_request = _comm.irecv(source = self.rank, tag = _TAG_PROGRESS)
		(has_update, progress) = self._progress_request.test()
		if has_update:
			self.progress = progress
			self._progress_request = None
		return has_update

class MPIJobManager:
	"""
	A class to manage the running of MPI jobs.
	"""
	def __init__(self):
		log_diagnostic(f"World size is {_comm.size}")
		self._workers = [MPIWorker(x) for x in range(1, _comm.size)]
		log_diagnostic(f"Using {len(self._workers)} workers...")
		self._lock = BoundedSemaphore(1)
		self._completed: list[Job] = []
		self.total_weight: int = 0

	def _any_idle(self) -> bool:
		"""
		Check if any workers are idle.
		"""
		return any([x for x in self._workers if x.is_idle()])

	def _get_idle(self) -> MPIWorker:
		"""
		Get an idle worker.
		"""
		for worker in self._workers:
			if worker.is_idle():
				return worker
		raise ValueError(f"No idle workers available")

	def _submit_job(self, job: Job, worker: MPIWorker):
		"""
		Submit the specified job to the specified worker.

		@param job: The job to be run.
		@param worker: The worker on which the job should be executed.
		"""
		with self._lock:
			if worker.job is not None:
				log_diagnostic(f"Worker {worker.rank} has completed a job")
				self._completed.append(worker.job)
			worker.job = job
		_mpi_master_send_job(worker.rank, job.task)

	def run(self, jobs: list[Job]):
		"""
		Run the specified jobs.

		@param jobs: Jobs to be run.
		"""
		self.total_weight = sum(j.weight for j in jobs)
		log_diagnostic(f"{len(jobs)} jobs to be run")
		for job in jobs:
			# If number of running jobs exceeds available number of CPUs,
			# wait for one job to finish.
			if not self._any_idle():
				log_diagnostic(f"No workers are currently available.")
				self._wait_until(self._any_idle)

			worker: MPIWorker = None
			with self._lock:
				worker = self._get_idle()
				worker._idle = False

			log_diagnostic(f"Job will be submitted to worker {worker.rank}")

			self._submit_job(job, worker)

		# Wait until job starts, to prevent race conditions.
		self._wait_until(lambda: all([worker.is_idle() for worker in self._workers]))

		# Signal to the workers that they may exit.
		for worker in self._workers:
			_comm.send(_CMD_DONE, worker.rank, _TAG_CMD)

	def _get_progress(self) -> float:
		"""
		Get overall progress as a float in range [0, 1].
		"""
		progress = 0
		for job in self._completed:
			progress += job.weight
		for worker in self._workers:
			if worker.job is not None:
				progress += worker.progress * worker.job.weight
		return progress / self.total_weight

	def _write_progress(self):
		"""
		Called after receiving a progress update from the specified worker.
		"""
		progress = self._get_progress()
		# log_diagnostic(f"Updating progress report with progress = {progress}")
		log_progress(progress)

	def _wait_until(self, condition: Callable[[None], bool]):
		"""
		Read progress messages from jobs until the given condition returns true.

		@param condition: The exit condition for this function.
		"""
		while True:
			if condition():
				return
			for worker in self._workers:
				if not worker.is_idle() and worker.has_progress_update():
					self._write_progress()
				if condition():
					return

def _mpi_worker_run():
	log_diagnostic(f"Worker is ready")
	while True:
		log_diagnostic(f"Sending idle notification")
		_comm.send(_comm.rank, _RANK_MASTER, _TAG_WORKER_IDLE)
		log_diagnostic(f"Waiting for command")
		cmd = _comm.recv(source = _RANK_MASTER, tag = _TAG_CMD)
		log_diagnostic(f"Received command from master")
		if cmd == _CMD_DONE:
			log_diagnostic(f"Received CMD_DONE: exiting")
			break
		elif cmd == _CMD_JOB:
			# Get the task to be executed.
			log_diagnostic(f"Received CMD_JOB: requesting job details")
			task: Task = _comm.recv(source = _RANK_MASTER, tag = _TAG_TASK)
			log_diagnostic(f"Received job from master. Job will now execute")
			try:
				task.exec(_mpi_worker_report_progress)
				log_diagnostic(f"Job completed successfully.")
			except BaseException as err:
				msg = str.join("", format_exception(err))
				log_error(f"Job ran with error:")
				log_error(msg)
				log_diagnostic(f"Sending error notification to master")
				_comm.send(msg, _RANK_MASTER, _TAG_ERROR)
		else:
			raise ValueError(f"Unknown command received from master: {cmd}")
	sys.exit(0)

def mpi_init():
	"""
	Initialise the MPI environment. This essentially acts as the main function
	if this is an MPI worker node, in which case control will NOT return from
	this function.
	"""
	comm = mpi4py.MPI.COMM_WORLD
	if comm.rank != _RANK_MASTER:
		_mpi_worker_run()
