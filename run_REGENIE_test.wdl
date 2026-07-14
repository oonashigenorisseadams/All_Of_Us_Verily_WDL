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
    set -uo pipefail

    echo "=== Attempting apt-get install libgomp1 (30s/60s timeouts) ==="
    timeout 30 apt-get update -qq
    UPDATE_RC=$?
    echo "apt-get update exit code: ${UPDATE_RC}"

    if [ "${UPDATE_RC}" -eq 0 ]; then
      timeout 60 apt-get install -qq -y --no-install-recommends libgomp1
      INSTALL_RC=$?
      echo "apt-get install exit code: ${INSTALL_RC}"
    else
      echo "Skipping install since apt-get update failed/timed out."
    fi

    echo "=== Checking for libgomp.so.1 ==="
    find / -name "libgomp.so*" 2>/dev/null || echo "libgomp.so not found anywhere"

    echo "=== regenie_bin ==="
    ls -lh ~{regenie_bin}
    chmod 750 ~{regenie_bin}

    echo "=== regenie version ==="
    ~{regenie_bin} --version || echo "no --version output"

    echo "=== regenie help ==="
    ~{regenie_bin} --help || echo "regenie --help failed, exit code $?"
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
