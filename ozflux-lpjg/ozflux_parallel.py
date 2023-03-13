from typing import Callable
import multiprocessing, os
import ozflux_logging

class _Job(multiprocessing.Process):
	"""
	Represents a parallel job.
	"""
	def __init__(self, id: int, func: Callable[[Callable[[float], None]], None]
	    	, progress_writer: multiprocessing.connection.Connection):
		"""
		Create a new _Job instance.

		@param id: Unique ID assigned to this job (used for progress reporting).
		@param func: The function to be called which performs work.
		@param progres_writer: Connection to the subprocess' stdout pipe.
		"""
		multiprocessing.Process.__init__(self)
		self.id = id
		self.func = func
		self.progress_writer = progress_writer

	def run(self):
		# Do main processing.
		self.func(lambda p: self.progress_writer.send( (p, self.id) ))

		# Close progress reporter pipe.
		self.progress_writer.close()

class JobManager:
	def __init__(self):
		# Global job weights. Access to this should only occur via _weights_lock.
		self._job_weights: list[int] = []

		# Total weight of all jobs. Access to this is controlled by _weights_lock.
		self._total_weight: int = 0

		# List of readable connections to subprocess' stdout pipes.
		self._readers: list[multiprocessing.connection.Connection] = []

		# List of progress for each job.
		self._overall_progress: list[float] = []

		# List of all submitted jobs (not running jobs!).
		self._jobs: list[_Job] = []

	def _progress_reporter(self, progress: float, id: int):
		"""
		This function is called by the wait() function when a progress report is
		received from one of the child processes.

		@param progress: Progress of this job in range [0, 1].
		@param id: ID of the job which has reported its progress.
		"""

		weight = self._job_weights[id] / self._total_weight
		job_progress = weight * progress

		# No need to check if in parallel mode, as the mutex is easily
		# obtained when running in serial mode.
		self._overall_progress[id] = job_progress
		aggregate_progress = sum(self._overall_progress)
		ozflux_logging.log_progress(aggregate_progress)

	def add_parallel_job(self, weight: int
		    , func: Callable[[Callable[[float], None]], None]):
		"""
		Add a new parallel job which uses the specified progress callback function.

		No progress reporting will occur until wait() is called.

		@param func: The function to be executed in parallel.
		@param weight: Absolute weight for this job. Higher value means progress in
					this job counts for proportionately more out of the aggregate
					progress of all jobs. The way that the weight is calculated
					should be consistent over all jobs, and it must be positive.
		"""
		# Get a job ID.
		job_id = len(self._jobs)

		# Update job weights. Note that the newly-added job weight should have
		# index job_id.
		self._job_weights.append(weight)
		self._total_weight += weight

		# Create a pipe for 1-way communication (progress reporting).
		reader, writer = multiprocessing.connection.Pipe(duplex = False)

		# Put the read pipe handle into the readers list.
		self._readers.append(reader)

		# Create an entry for this job in the progress list.
		self._overall_progress.append(0)

		# Start the process.
		job = _Job(job_id, func, writer)
		job.start()

		# Close the writable end of the pipe now, to be sure that p is the only
		# process which owns a handle for it. This ensures that when p closes its
		# handle for the writable end, wait() will promptly report the readable end
		# as being ready.
		writer.close()

		# Store the process handle for later use.
		self._jobs.append(job)

	def wait(self):
		"""
		Wait for all jobs to finish running. Note that no progress messages will be
		written until this is called.
		"""
		while self._readers:
			for reader in self._readers:
				try:
					msg = None
					if reader.poll():
						msg = reader.recv()
						(progress, process_index) = msg
						self._progress_reporter(progress, process_index)
				except EOFError:
					self._readers.remove(reader)

		# Wait for processes to exit (actually, they should have all finished by
		# the time we get to here).
		for process in self._jobs:
			process.join()
