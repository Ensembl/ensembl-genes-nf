Pipeline-owned orchestration scripts belong here when they are not directly
part of a tool process. Small deterministic data transformations live in
`bin/` so they can be unit tested and invoked reproducibly from a task.
