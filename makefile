compila : main_iplib.c bmp.o ip_lib.o
	gcc main_iplib.c -o main_iplib ip_lib.o bmp.o -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

bmp.o : bmp.c bmp.h
	gcc -c bmp.c -o bmp.o -Wall -lm
  
ip_lib.o : ip_lib.c ip_lib.h
	gcc -c ip_lib.c -o ip_lib.o -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

clean:	
		rm ip_lib.o bmp.o main_iplib
