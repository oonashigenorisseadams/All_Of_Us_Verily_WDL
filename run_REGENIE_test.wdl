version 1.0

task regenie_help_test {
  input {
    File regenie_bin
    File bedfile
    File bimfile
    File famfile
    File phenocovar_file
    Int cpu
    Int mem
    Int disk_gb
  }

  command <<<
    set -euo pipefail

    echo "=== Staged files ==="
    find . -maxdepth 6 -type f

    echo "=== Disk usage ==="
    df -h .

    echo "=== Memory ==="
    free -h

    echo "=== bedfile ==="
    ls -lh ~{bedfile}

    echo "=== bimfile ==="
    ls -lh ~{bimfile}

    echo "=== famfile ==="
    ls -lh ~{famfile}

    echo "=== phenocovar_file ==="
    ls -lh ~{phenocovar_file}

    echo "=== regenie_bin ==="
    ls -lh ~{regenie_bin}
    chmod 750 ~{regenie_bin}

    echo "=== regenie version ==="
    ~{regenie_bin} --version || echo "no --version output"

    echo "=== regenie help ==="
    ~{regenie_bin} --help
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

workflow test_regenie_localization {
  input {
    File regenie_bin
    File bedfile
    File bimfile
    File famfile
    File phenocovar_file
    Int cpu = 4
    Int mem = 25
    Int disk_gb = 300
  }

  call regenie_help_test {
    input:
      regenie_bin = regenie_bin,
      bedfile = bedfile,
      bimfile = bimfile,
      famfile = famfile,
      phenocovar_file = phenocovar_file,
      cpu = cpu,
      mem = mem,
      disk_gb = disk_gb
  }

  output {
    File stdout_log = regenie_help_test.stdout_log
    File stderr_log = regenie_help_test.stderr_log
  }
}
