# CFLAGS = -O4 -DDEBUG

default: app

app: main.c haralick_imp.h haralick_imp.c stb_image.h misc.c misc.h
	smpicc $(CFLAGS) main.c haralick_imp.c haralick_imp.h stb_image.h misc.c misc.h -o app

run:
	smpirun -np 4 -hostfile ./config/hostfile_64.txt -platform ./config/cluster_crossbar_64.xml ./app ./data/random.png 0
	smpirun -np 1 -hostfile ./config/hostfile_64.txt -platform ./config/cluster_crossbar_64.xml ./app ./data/random.png 0 -i sequential

clean:
	/bin/rm -f */*.o *.o smpitmp* app *.out *.err