version 1.0

task sum {
  input {
    Array[String] nums
  }
  
  command <<<
  printf ~{sep(" ", nums)} | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}'
  >>>
  
  runtime {
    docker: 'ubuntu:latest'
  }
  output {
    Int total = read_int(stdout())
  }
}