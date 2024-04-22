from typing import Callable
from multiprocessing import BoundedSemaphore
from ozflux_logging import *
import concurrent.futures, mpi4py, mpi4py.futures
from ozflux_parallel import Task, Job

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

# Tag used for transmission of number of arguments from master to worker nodes.
_TAG_NARG = 6

# Tag used for transmission of job arguments from master to worker nodes.
_TAG_ARG = 7

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
	_comm.send(task.exec, worker, _TAG_TASK)

class MPIWorker:
	"""
	Stores metadata and handles communication with an MPI worker.
	"""
	def __init__(self, rank: int):
		self.rank = rank
		self._idle = True
		self._progress: float = 0.0
		self._idle_request: mpi4py.MPI.Request = None
		self._progress_request: mpi4py.MPI.Request = None

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
			self._idle_request = _comm.irecv(source = self.rank, tag = _TAG_WORKER_IDLE)

		# If the request for a WORKER_IDLE message has completed, the worker is
		# now idle.
		if self._idle_request.test():
			self._idle_request = None
			self._idle = True

		return self._idle

	def get_progress(self) -> float:
		"""
		Get the progress of this job as a fraction in range [0, 1].
		"""
		if self._idle:
			return 1.0

		if self._progress_request is None:
			self._progress_request = _comm.irecv(source = self.rank, tag = _TAG_PROGRESS)

		if self._progress_request.test():
			self._progress = self._progress_request.wait()
			self._progress_request = None

		return self._progress

class MPIJobManager:
	"""
	A class to manage the running of MPI jobs.
	"""
	def __init__(self):
		self._workers = [MPIWorker(x) for x in range(1, _comm.size)]
		self._lock = BoundedSemaphore(1)

	def _get_num_busy(self) -> int:
		"""
		Get the number of busy workers.
		"""
		return len([x for x in self._workers if not x._idle])

	def _any_idle(self) -> bool:
		"""
		Check if any workers are idle.
		"""
		return any([x for x in self._workers if x._idle])

	def _get_idle(self) -> MPIWorker:
		"""
		Get an idle worker.
		"""
		for worker in self._workers:
			if worker._idle:
				return worker
		raise ValueError(f"No idle workers available")

	def run(self, jobs: list[Job]):
		"""
		Run the specified jobs.

		@param jobs: Jobs to be run.
		"""
		for job in jobs:
			# If number of running jobs exceeds available number of CPUs,
			# wait for one job to finish.
			if not self._any_idle():
				self._wait_until(self._any_idle)

			worker: MPIWorker = None
			with self._lock:
				worker = self._get_idle()
				worker._idle = False

			_mpi_master_send_job(worker, job.task)

		# Wait until job starts, to prevent race conditions.
		self._wait_until(lambda: len(self._idle) == _comm.size - 1)

	def _wait_until(self, condition: Callable[[float], None]):
		"""
		Read progress messages from jobs until the given condition returns true.

		@param condition: The exit condition for this function.
		"""
		for worker in self._workers:
			if not worker.is_idle():
				if 

		while workers:
			for job in jobs:
				try:
					if condition():
						ozflux_logging.log_diagnostic(f"All jobs appear to have finished running")
						return
					if not job.is_running():
						ozflux_logging.log_diagnostic(f"Job {job.get_id()} appears to have finished running")
						self._progress_reporter(1, job.get_id())
						jobs.remove(job)
						continue
					if job.has_progress_update():
						ozflux_logging.log_diagnostic(f"Job {job.get_id()} has a progress update")
						self._progress_reporter(job.get_progress(), job.get_id())
				except EOFError:
					ozflux_logging.log_diagnostic(f"EOFError from job {job.get_id()}")
					jobs.remove(job)
		ozflux_logging.log_diagnostic(f"All jobs removed from queue")

def _mpi_worker_run():
	while True:
		_comm.send(_comm.rank, _RANK_MASTER, _TAG_WORKER_IDLE)
		cmd = _comm.recv(source = _RANK_MASTER, tag = _TAG_CMD)
		if cmd == _CMD_DONE:
			break
		elif cmd == _CMD_JOB:
			# Get the task to be executed.
			task: Task = _comm.recv(source = _RANK_MASTER, tag = _TAG_TASK)
			task.exec(_mpi_worker_report_progress)
		else:
			raise ValueError(f"Unknown command received from master: {cmd}")
	sys.exit(0)

def mpi_init():
	comm = mpi4py.MPI.COMM_WORLD
	if comm.rank != _RANK_MASTER:
		_mpi_worker_run()
