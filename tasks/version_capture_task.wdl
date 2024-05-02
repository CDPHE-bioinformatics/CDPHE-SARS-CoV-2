version 1.0

struct VersionInfo {
  String software
  String docker
  String version
}

# workaround cromwell bug with read_json of Array
# https://github.com/openwdl/wdl/issues/409
struct VersionInfoArray {
  Array[VersionInfo] versions
}

task workflow_version_capture {
  input {
    String? timezone
  }
  meta {
    description: "capture version release"
  }
  command <<<

    Workflow_Version="v2-3-0"

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

task task_version_capture {
  input {
    Array[VersionInfo] version_array
    String workflow_name
    String workflow_version
    String project_name
    String analysis_date
    File version_capture_py
  }

  VersionInfoArray versions = { "versions": version_array }

  command <<<
    python ~{version_capture_py} \
      --versions_json ~{write_json(versions)} \
      --workflow_name ~{workflow_name} \
      --workflow_version ~{workflow_version} \
      --project_name ~{project_name} \
      --analysis_date ~{analysis_date}
  >>>

  output {
    File version_capture_file = "version_capture_${workflow_name}_${project_name}_${workflow_version}.csv"
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "mchether/py3-bio:v4"
    disks: "local-disk 10 SSD"
  }
}
