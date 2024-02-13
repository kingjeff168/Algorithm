solution: main.o
	gcc -o particles main.o -lm

%.o: %.c
	gcc -c -o $@ $<

clean:
	rm -f *.o particles