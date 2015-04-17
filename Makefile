IC       = icc
GC       = g++
CFLAS   = -tpp2
OPTFLAGS = -O3 -static -fno-alias
PROFFLAGS= -pg -O0
PROFEXEC = Model_prof
EXEC     = Model
SCENFLAGS= -DDEBUG_ON -DW_SYSDATA -DCLIM_SCEN #-DS_FIRE_FROM_FILE

all:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA  #-DOPTIM_AUSTRALIA #-DS_FIRE_FROM_FILE
# 	/usr/bin/time ./Model 100 150 31.266  -25.166 0 0 0 test  # pkop
# 	/usr/bin/time ./Model 100 150 31.266  -25.166 1 0 0 test  # pkop
# 	/usr/bin/time ./Model 100 150 31.77   -24.39  1 0 0 test  # satara
# 	/usr/bin/time ./Model 100 500 31.77   -24.39  0 0 0 test  # satara
	/usr/bin/time ./Model 100 200 31.58  -24.99  0 0 0 test  # skukuza
# 	/usr/bin/time ./Model 100 300 18.5  -34  1 0 0 test  # fynbos

# 	/usr/bin/time ./Model 100 200 31.58   -24.99  1 0 0 test  # skukuza
# 	/usr/bin/time ./Model 100 200 31.58   -24.99  0 0 0 test  # skukuza



opti:
	${GC} -o Model Model.cpp -DDEBUG_ON -DOPTIM_GLOBALS
	./Model 100 200 31.58  -24.99  0 0 0 test 0.001 0.001 0.0001 0.0005 0.001 0.01 0.001 0.001 0.001 0.01 0.1 0.5 0.5 0.5 0.5 0.5 0.5 0.9 0.9 0.2 0.2 0.5 0.9 0.2 0.5 0.9 0.2 0.2 0.2 0.2 0.5 0.5 0.9 0.9


sav0:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 31.58 -24.99 1 0 3 test
# 	./Model 1000 339 31.58 -24.99 1 0 3 test

les0:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 28 -30 1 0 3 test

les1:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 30 -27 1 0 3 test

fyn0:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 19 -34   1 0 3 test
# 	./Model 1000 339 19 -34   1 0 3 test


fyn1:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 19.583 -34.583 1 0 3 test
# 	./Model 1000 339 19     -34   1 0 3 test

fyn2:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 18.417 -33.917 1 0 3 test

ecp:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 30 -30 1 0 3 test
# 	./Model 1000 339 30 -30 1 0 3 test


fvar:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 1000 339 19 -34.0001  1 0 3 test

rf:
	${GC} -o Model Model.cpp ${SCENFLAGS}
	./Model 100 339 20 0         1 0 3 test


test:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA -DS_NOREPROD -DS_NOGRASS
	./Model 3 60 31.58  -24.99   0 0 3 test | grep LLL > sav0 &
	./Model 3 60 28     -30      0 0 3 test | grep LLL > les0 &
	./Model 3 60 19     -34      0 0 3 test | grep LLL > fyn0 &
	./Model 3 60 19.583 -34.583  0 0 3 test | grep LLL > fyn1 &
	./Model 3 60 18.417 -33.917  0 0 3 test | grep LLL > fyn2 &
	./Model 3 60 19     -34.0001 0 0 3 test | grep LLL > fynf &

csens:
# 	${GC} -o Model Model.cpp -DDEBUG_ON -DS_CLIM_SENSITIVITY -DW_SYSDATA -DS_CLIM_NO_C4
	${GC} -o Model Model.cpp -DDEBUG_ON -DS_CLIM_SENSITIVITY -DW_SYSDATA 
	/usr/bin/time ./Model 100 300 31.58   -24.99  1 0 0 test 387 0 1 0 # skukuza

cc:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA -DCLIM_SCEN #-DOPTIM_GLOBALS #-DS_FIRE_FROM_FILE
# 	/usr/bin/time ./Model 100 339 31.266  -25.166 0 0 3 test  # pkop
# 	/usr/bin/time ./Model 100 339 31.266  -25.166 1 0 3 test  # pkop
# 	/usr/bin/time ./Model 100 339 31.77   -24.39  1 0 3 test  # satara
# 	/usr/bin/time ./Model 100 339 31.77   -24.39  0 0 3 test  # satara
	/usr/bin/time ./Model 100 339 31.58   -24.99  1 0 3 test  # skukuza
# 	/usr/bin/time ./Model 100 339 31.58   -24.99  0 0 3 test  # skukuza

aus:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA -DS_AUSTRALIA #-DS_FIRE_FROM_FILE
# 	/usr/bin/time ./Model 100 150 131.15  -12.49  0 0 0 test
	/usr/bin/time ./Model 100 150 131.15  -12.49  1 0 0 test

bra:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA #-DS_FIRE_FROM_FILE
# 	/usr/bin/time ./Model 100 250 -60.017 -3.133  0 0 0 test
# 	/usr/bin/time ./Model 100 250 -60.017 -3.133  1 0 0 test
# 	/usr/bin/time ./Model 100 250 -54.969 -3.018  1 0 0 test
	/usr/bin/time ./Model 100 250 -54.969 -3.018  0 0 0 test

prof:
	${GC} -o ${PROFEXEC}  ${PROFFLAGS} Model.cpp

#prof:
#	${IC} -o ${PROFEXEC} ${CFLAGS} ${OPTFLAGS} ${PROFFLAGS} Model.cpp

opt: 
	${IC} -o ${EXEC} ${CFLAGS} ${OPTFLAGS} Model.cpp

dbg: 
	${IC} -o ${EXEC} ${CFLAGS} ${PROFFLAGS} -DDEBUG Model.cpp

g++: 
	${GC} -o ${EXEC} ${CFLAGS} ${OPTFLAGS} Model.cpp

profile:
	icc  -O3 -pg -o Model_profile Model.cpp
	./Model_profile 1500  100    -18.18    31.47        0   0   0     test

IDataMaps:
	icc  -O3  -fno-fnalias -fno-exceptions -static -o Model_IDataMaps Model.cpp  -DW_IDATAMAPS


ele:
	${GC} -o Model Model.cpp -DDEBUG_ON -DW_SYSDATA -DS_ELEPHANTS #-DOPTIM_AUSTRALIA #-DS_FIRE_FROM_FILE
# 	/usr/bin/time ./Model 100 150 31.266  -25.166 0 0 0 test  # pkop
	/usr/bin/time ./Model 100 500 31.266  -25.166 1 0 0 test 5 5 0.1 0.1 0.5 0.1 12000 0.5 150 0.1 # pkop
# 	/usr/bin/time ./Model 100 150 31.77   -24.39  1 0 0 test  # satara
# 	/usr/bin/time ./Model 100 500 31.77   -24.39  0 0 0 test  # satara
# 	/usr/bin/time ./Model 100 500 31.58   -24.99  1 0 0 test  # skukuza
# 	/usr/bin/time ./Model 100 200 31.58   -24.99  0 0 0 test  # skukuza



run:
	date
	/usr/bin/time ./Model 100 600 -24.6 28.7 0 1 0
	date

# # profile:
# # 	date
# # 	profile.pl -c5 ./Model 100 600 -24.6 28.7 0 1 0
# # 	date

clean:
	rm -rf ${EXEC}

