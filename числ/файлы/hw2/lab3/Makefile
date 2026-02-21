lab3: wip.c
	gcc -std=c23 \
	-Wall -Wextra -Wpedantic \
	-Wconversion -Wshadow -Wstrict-prototypes -Wfloat-equal -fno-common \
	-g3 -Og \
	-fsanitize=address,undefined -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fstack-protector-all \
	wip.c -o lab3 -lm -llapacke -lopenblas
	./lab3
