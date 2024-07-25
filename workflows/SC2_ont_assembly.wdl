version 1.0

import "../tasks/repo_version.wdl" as repo_version_task


workflow SC2_ont_assembly {
    call repo_version_task.print_version as repo_version {}
    output {
        String version = repo_version.repo_version_output
    }
}