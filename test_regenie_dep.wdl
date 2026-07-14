version 1.0

task check_ubuntu {
  input {
    File regenie_bin
    Int cpu
    Int mem
    Int disk_gb
  }

  command <<<
    set -uo pipefail

    echo "=== ldd on regenie_bin under ubuntu:22.04 ==="
    ldd ~{regenie_bin} 2>&1 || echo "(ldd returned non-zero, see output above)"

    echo "=== Attempting to run regenie anyway ==="
    chmod 750 ~{regenie_bin}
    ~{regenie_bin} --version 2>&1 || echo "no --version output, exit $?"
    ~{regenie_bin} --help 2>&1 || echo "regenie --help failed, exit code $?"
  >>>

  output {
    File stdout_log = stdout()
    File stderr_log = stderr()
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task check_gcc {
  input {
    File regenie_bin
    Int cpu
    Int mem
    Int disk_gb
  }

  command <<<
    set -uo pipefail

    echo "=== ldd on regenie_bin under gcc:12 ==="
    ldd ~{regenie_bin} 2>&1 || echo "(ldd returned non-zero, see output above)"

    echo "=== Attempting to run regenie ==="
    chmod 750 ~{regenie_bin}
    ~{regenie_bin} --version 2>&1 || echo "no --version output, exit $?"
    ~{regenie_bin} --help 2>&1 || echo "regenie --help failed, exit code $?"
  >>>

  output {
    File stdout_log = stdout()
    File stderr_log = stderr()
  }

  runtime {
    docker: "gcc:12"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

workflow test_regenie_deps {
  input {
    File regenie_bin
    Int cpu = 4
    Int mem = 25
    Int disk_gb = 300
  }

  call check_ubuntu {
    input:
      regenie_bin = regenie_bin,
      cpu = cpu,
      mem = mem,
      disk_gb = disk_gb
  }

  call check_gcc {
    input:
      regenie_bin = regenie_bin,
      cpu = cpu,
      mem = mem,
      disk_gb = disk_gb
  }

  output {
    File ubuntu_stdout = check_ubuntu.stdout_log
    File ubuntu_stderr = check_ubuntu.stderr_log
    File gcc_stdout = check_gcc.stdout_log
    File gcc_stderr = check_gcc.stderr_log
  }
}
