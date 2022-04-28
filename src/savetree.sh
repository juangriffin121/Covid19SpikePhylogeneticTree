#!/bin/bash
echo "_________"
echo "Make Tree"
echo "_________"
echo "press 1 if you will enter your sequence as a text file, 2 to enter your sequences in the terminal or 0 to exit the program"
from_txt(){
	./make_tree.py $1 $2
}

from_input(){
        echo "insert your sequence in FASTA format enter a "." at the end of your last sequence"
	read -d "." sequences
		echo $sequences > temp
		from_txt temp arbol
		rm temp
}

insert_txt(){
	echo "insert the path to your file"
	read path
		from_txt $path arbol
}

read -n 1 -s choice;
	case $choice in
		2)from_input;;
		1)insert_txt;;
		0) exit;;
	esac
