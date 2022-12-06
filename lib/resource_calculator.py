def compute_octopus_mem_mb(wildcards, input, threads, attempt) -> int:
    """
    Compute scaling RAM requirements for octopus_run_task with each resubmission
    to the compute controller
    """
    octopus_scaling_ram = 16000
    return int(threads * octopus_scaling_ram * (1 + (attempt - 1) * 0.1))
