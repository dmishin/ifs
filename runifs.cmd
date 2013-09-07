@echo off

echo compile
g++ -Wall -O3 ifs-rand.cpp -o ifs-rand-sm.exe || goto :error

echo run
ifs-rand-sm.exe                     || goto :error 

echo show result
start test-small.pgm                   
goto :EOF

:error
echo Error occured.
     