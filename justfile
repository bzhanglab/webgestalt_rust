test: test-exp test-lib
  @echo "Done Testing"
test-exp:
    cargo run -F "experimental" -- integrated --help

[working-directory('webgestalt_lib')]
test-lib:
    cargo test
