final: main.c mesh.c aux.c
	gcc -Wall main.c mesh.c aux.c -lm -o shockTube
