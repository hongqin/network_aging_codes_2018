if [ -f Rplots.ps ]
then
 rm Rplots.ps
fi
R --no-save < $1
if [ -f Rplots.ps ]  
then
 gs Rplots.ps &
fi
exit 0
