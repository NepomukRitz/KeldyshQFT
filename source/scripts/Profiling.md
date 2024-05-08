## Profiling with gprof
1. compile with flag -pg (in CMakeLists.txt)
2. running executable  gives gmon.out
3. run gprof
   1. grof <executable> gmon.out > analysis.text
   2. After downloading grof2dot:
   gprof <executable> gmon.out | gprof2dot-master/gprof2dot.py > graph.dot