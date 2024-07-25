version 1.0

task print_version {
    command <<<
        git describe --tags | tee repo_version_str
    >>>
    output {
        String repo_version_output = read_string("repo_version_str")
    }
}