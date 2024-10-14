# R test SLURM arrary task id
print(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))