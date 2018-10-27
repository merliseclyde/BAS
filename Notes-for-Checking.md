# To check with valgrind on fedora
#
R -d "valgrind --tool=memcheck --leak-check=full --log-fd=1 --log-file=BAS.Rcheck/BAS-valgrind.txt" --vanilla < BAS.Rcheck/BAS-Ex.R
