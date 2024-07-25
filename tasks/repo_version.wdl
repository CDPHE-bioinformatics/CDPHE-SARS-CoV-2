version 1.0

workflow version_testing {
    call print_version {}
}
task print_version {
    command <<<
        git describe --tags | tee repo_version
    >>>
    output {
        String repo_version_output = read_string("repo_version")
    }
}