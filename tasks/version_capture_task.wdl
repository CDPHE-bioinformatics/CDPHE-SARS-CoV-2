version 1.0

task workflow_version_capture {
  input {
    String? timezone
  }
  meta {
    description: "capture version release"
  }
  command <<<
    Workflow_Version="v2-1-2"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$Workflow_Version" > WORKFLOW_VERSION
  >>>
  output {
    String analysis_date = read_string("TODAY")
    String workflow_version = read_string("WORKFLOW_VERSION")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu:jammy"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}