#! /bin/bash
#SBATCH --job-name="Hello"
#SBATCH --gpus=0
#SBATCH --nodes=1                      # Количество используемых узлов
#SBATCH --ntasks=1                     # Количество MPI процессов
#SBATCH --cpus-per-task=44             # Количество ядер CPU, выделенных для одного процесса
#SBATCH --cores-per-socket=22          # Количество ядер на одном сокете узла
#SBATCH --time=0-0:20
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_b"

# Executable
srun ./hello
