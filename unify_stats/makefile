##################################################################
#                        Makefile for f90                         #
###################################################################

########################     compiler     #########################

F90 = ifort
F77 = ifort
#F90 = gfortran 
#F77 = gfortran 

########################  compiler flags  #########################

F90FLAGS= -c -O3
#F90FLAGS = -c -warn -CB -debug extended -traceback
F77FLAGS= -c -O3
LFLAGS =
#PREP = scalasca -inst
#PREP = skin

########################  objects alpha   #########################

INIT    = .
SRCDIR  = $(INIT)
OBJ     = $(INIT)
OBJDIR  = $(OBJ)
CALCDIR = $(INIT)

OBJECTS = $(OBJ)/declaration_unify.o\
          $(OBJ)/unify_stats.o
          
########################      build       #########################

unify_stats :printmsgA $(OBJECTS)
	@echo
	@echo Linking... 
	$(PREP) $(F90) -o $@ $(OBJECTS) $(LFLAGS)
	@echo 
	@echo unify_stats built, congratulations.
	@echo 

########################     compile      #########################

$(OBJDIR)/declaration_unify.o : $(SRCDIR)/declaration_unify.f90 $(SRCDIR)/makefile
	@echo compiling declaration_unify.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/declaration_unify.f90

$(OBJDIR)/unify_stats.o : $(SRCDIR)/unify_stats.f90  $(SRCDIR)/makefile
	@echo compiling unify_stats.f90
	@cd $(OBJDIR); $(PREP) $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/unify_stats.f90

########################      message     #########################

printmsgA :
	@echo
	@echo Building ...
	@echo F77 Compiler flags : $(F77FLAGS)
	@echo F90 Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)
	@echo Prepend  flags : $(PREP)

########################      clean       #########################

clean:
	@find . \( -name '*.o' \) -exec rm {} \;

########################   end of file    #########################
