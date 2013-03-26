
CC = g++ -o3 

#INCLUDES = -I/home/risi/ode/include 
# -I/usr/include/glib-1.2 3-I/usr/include/g++-2 -I/usr/lib/gtkmm/include -I/usr/lib/sigc++/include -I/usr/lib/glib/include -I/usr/include/gtk-1.2

#CFLAGS = -g -Wall -Wno-return-type $(INCLUDES) -DSWIG_GLOBAL
#CFLAGS = -g -Wall -Werror
CFLAGS = -g -Wall 

biped: biped.o neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o organism.o species.o population.o evolvable_substrate.o noveltyset.o
	$(CC) $(CFLAGS) -I/home/risi/ode/include -L/home/risi/ode/lib biped.o neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o noveltyset.o organism.o evolvable_substrate.o species.o population.o -lode -o biped
#	$(CC) $(CFLAGS) $(LIBS) networks.o genetics.o -o neat `gtkmm-config --cflags --libs`

########################
biped.o: biped.cpp neat.h
	$(CC) $(CFLAGS) -I/home/risi/ode/include -c biped.cpp -o biped.o

neat.o: neat.cpp neat.h
	  $(CC) $(CFLAGS) -c neat.cpp -o neat.o

network.o: network.cpp network.h neat.h neat.o  
	$(CC) $(CFLAGS) -c network.cpp -o network.o

nnode.o: nnode.cpp nnode.h    
	$(CC) $(CFLAGS) -c nnode.cpp -o nnode.o

link.o: link.cpp link.h
	  $(CC) $(CFLAGS) -c link.cpp -o link.o

trait.o: trait.cpp trait.h
	  $(CC) $(CFLAGS) -c trait.cpp -o trait.o

gene.o: gene.cpp gene.h
	  $(CC) $(CFLAGS) -c gene.cpp -o gene.o

genome.o: genome.cpp genome.h evolvable_substrate.h
	  $(CC) $(CFLAGS) -c genome.cpp -o genome.o

evolvable_substrate.o: evolvable_substrate.cpp evolvable_substrate.h
	$(CC) $(CFLAGS) -c evolvable_substrate.cpp -o evolvable_substrate.o

innovation.o: innovation.cpp innovation.h
	  $(CC) $(CFLAGS) -c innovation.cpp -o innovation.o

organism.o: organism.cpp organism.h noveltyset.h
	$(CC) $(CFLAGS) -c organism.cpp -o organism.o

species.o: species.cpp species.h organism.h
	  $(CC) $(CFLAGS) -c species.cpp -o species.o

population.o: population.cpp population.h organism.h neat.h neat.o
	   $(CC) $(CFLAGS) -c population.cpp -o population.o

noveltyset.o: noveltyset.cpp noveltyset.h
	  $(CC) $(CFLAGS) -c noveltyset.cpp -o noveltyset.o

########################

clean:
	rm -f neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o organism.o species.o population.o experiments.o neatmain.o neat
