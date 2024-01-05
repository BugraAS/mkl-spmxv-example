# Building

you need to have mkl installed and initialized. Alongside open mp.  
Also you need to create a directory named `build` for compilation.  
  
Configure project.  
```bash
cmake -Ssrc -Bbuild
```

Compile the project
```bash
cmake --build build --config Release --target all
```

Run
```bash
./build/intel-spmxv <Filepath-to-matrix> <number-of-multiplications>
```
