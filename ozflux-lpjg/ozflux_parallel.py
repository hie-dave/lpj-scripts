from typing import Callable
import multiprocessing
import ozflux_logging, ozflux_mpi
import concurrent.futures
from ozflux_common import Task, Job

class LocalJob:
	"""
	Represents a job that is runnin the local execution environment.
	"""
	def __init__(self, job: Job):
		"""
		Create a new LocalJob instance.

		@param id: Job ID.
		@param weight: Weight of the job.
		@param task: The task to be executed.
		"""
		self._job = job
		self._start = 0.0
		self._total_weight = 0.0

	def _progress_update(self, progress: float):
		"""
		Write a progress update to stdout.
		"""
		# Update this job's progress.
		self._progress = progress

		# Calculate overall progress (for all jobs).
		overall = self._start + progress * self._job.weight / self._total_weight

		# Write progress update.
		ozflux_logging.log_progress(overall)

	def get_id(self) -> int:
		"""
		Get the job's ID.
		"""
		return self._job.id

	def set_weight(self, start: int, total_weight: int):
		"""
		Initialise the job weighting. This is required in order for overall
		progress to be calculated correctly.

		@param start: The weight of all jobs completed before this one.
		@param total_weight: The total weight of all jobs.
		"""
		self._start = start
		self._total_weight = total_weight

	def run(self) -> bool:
		"""
		Run the job. This will execute the job in the current thread and block
		until it is finished. Any exceptions will be propagated out of this
		function. Progress reports will be written directly to stdout as
		start + progress * weight / total_weight.
		"""
		self._job.task.exec(self._progress_update)

class ProcessJob(multiprocessing.Process):
	"""
	Represents a parallel job which will be run as an independent process.
	"""
	def __init__(self, job: Job):
		"""
		Create a new ProcessJob instance.

		@param job: The job to be run.
		"""
		multiprocessing.Process.__init__(self)
		self._job = job

		# Create a pipe for 1-way communication (progress reporting).
		reader, writer = multiprocessing.connection.Pipe(duplex = False)

		self._reader = reader
		self.writer = writer

	def get_id(self) -> int:
		"""
		Get the job's ID.
		"""
		return self._job.id

	def run(self):
		"""
		Run the job (as a separate process).
		"""
		# Progress callback function.
		pcb = lambda p: self.writer.send( (p, self._job.id) )
		self._job.task.exec(pcb)

		# Close progress reporter pipe.
		self.writer.close()

	def is_running(self) -> bool:
		"""
		Check if this job is currently running.
		"""
		# Check if the process is running.
		return self.is_alive()

	def has_progress_update(self) -> bool:
		"""
		Check if this job has a progress update.
		"""
		return self._reader.poll()

	def get_progress(self) -> float:
		"""
		Get the current progress of this job. This will block until a progress
		message is written from the worker process, so it must be called only
		after has_progress_update() returns true.
		"""
		msg = self._reader.recv()
		(progress, _) = msg
		return progress

class JobManager:
	def __init__(self):
		"""
		Create a new JobManager.
		"""
		# # True iff jobs are allowed to run in parallel. False otherwise.
		# self.allow_parallel = allow_parallel

		# Total weight of all jobs. Access to this is controlled by _weights_lock.
		self._total_weight: int = 0

		# List of all submitted jobs (not running jobs!).
		self._jobs: list[Job] = []

		self._lock = multiprocessing.BoundedSemaphore(1)

	def _progress_reporter(self, progress: float, id: int):
		"""
		This function is called by the wait() function when a progress report is
		received from one of the child processes.

		@param progress: Progress of this job in range [0, 1].
		@param id: ID of the job which has reported its progress.
		"""

		job = self._jobs[id]

		weight = job.weight / self._total_weight

		job.progress = weight * progress

		# No need to check if in parallel mode, as the mutex is easily
		# obtained when running in serial mode.
		aggregate_progress = sum([j.progress for j in self._jobs])
		ozflux_logging.log_progress(aggregate_progress)

	def add_job(self, task: Task, weight: int = 1):
		"""
		Register a job with the job manager. A job can be any function which
		reports its progress.

		No progress reporting will occur until wait() is called.

		@param task: The job's execution function. This is any function which
					 reports progress via a callable argument.
		@param weight: Absolute weight for this job. Higher value means progress
					in this job counts for proportionately more out of the
					aggregate progress of all jobs. The way that the weight is
					calculated should be consistent over all jobs, and it must
					be positive.
		"""
		with self._lock:
			# Get a job ID.
			job_id = len(self._jobs)

			# Update job weights. Note that the newly-added job weight should have
			# index job_id.
			self._total_weight += weight

			# Start the process.
			job = Job(job_id, weight, task)

			# Store the process handle for later use.
			self._jobs.append(job)

	def run_parallel(self, max_para: int = multiprocessing.cpu_count()):
		"""
		Run all jobs in parallel and wait for them to finish.
		@param max_para: Maximum number of jobs to run in parallel. Defaults to logical CPU count.
		"""
		with self._lock:
			processes: list[ProcessJob] = []
			get_num_jobs = lambda: len([p for p in processes if p.is_running()])

			for job in self._jobs:
				process = ProcessJob(job)
				processes.append(process)

				# If number of running jobs exceeds available number of CPUs,
				# wait for one job to finish.
				num_running = self._get_num_running_jobs()
				if num_running >= max_para:
					self._wait_until(lambda: get_num_jobs() < max_para)

				process.start()

				# Close the writable end of the pipe now, to be sure that p is
				# the only process which owns a handle for it. This ensures that
				# when p closes its handle for the writable end, wait() will
				# promptly report the readable end as being ready.
				process.writer.close()

				# Wait until job starts, to prevent race conditions.
				self._wait_until(lambda: process.is_alive())

			# Wait until all jobs are finished.
			self._wait_until()
			for process in processes:
				process.join()

	def run_single_threaded(self):
		"""
		Run all jobs one at a time, in the current thread, and wait for them to
		finish.
		"""
		with self._lock:
			cum_weight = 0
			for task in self._jobs:
				# Create a local job to manage the execution of this task.
				job = LocalJob(task)
				job.set_weight(cum_weight, self._total_weight)

				# Run the job.
				job.run()
				cum_weight += task.weight

	def run_mpi(self):
		"""
		Run jobs in parallel, using MPI.
		"""
		# if mpi4py.MPI.COMM_WORLD.size < 2:
		# 	raise ValueError(f"Cannot run with MPI: world size is <2")
		job_manager = ozflux_mpi.MPIJobManager()
		job_manager.run(self._jobs)

	def _wait_until(self, jobs, condition: Callable[[], bool] = lambda: False):
		"""
		Wait for all parallel jobs to finish running, or the given condition
		returns true. Progress updates will be written while waiting for jobs to
		finish running.

		@param condition: Optional function which returns a bool. If this
		function returns true, waiting will cease. If no condition is given,
		this function will wait until all jobs have finished.
		"""
		jobs = [j for j in jobs if j.is_running()]
		while jobs:
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
