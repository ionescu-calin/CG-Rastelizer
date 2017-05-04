FILE=skeleton

########
#   Directories
S_DIR=Source
B_DIR=Build

S_DIR_B=Source/Source_Basic
S_DIR_C=Source/Source_Clipping
S_DIR_F=Source/Source_Final
########
#   Output
EXEC=$(B_DIR)/$(FILE)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -fopenmp -ggdb -g3 
LN_OPTS=
CC=g++ -fopenmp -O3 

########
#       SDL options
SDL_CFLAGS := $(shell sdl-config --cflags)
GLM_CFLAGS := -I$(GLMDIR)
SDL_LDFLAGS := $(shell sdl-config --libs)

########
#   This is the default action
all:Build


########
#   Object list
#
OBJ = $(B_DIR)/$(FILE).o


########
#   Objects
$(B_DIR)/$(FILE).o : $(S_DIR)/$(FILE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)


########
#   Main build rule     
Build : $(OBJ) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

basic :  $(S_DIR_B)/$(FILE).cpp $(S_DIR_B)/SDLauxiliary.h $(S_DIR_B)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_B)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

clipping : $(S_DIR_C)/$(FILE).cpp $(S_DIR_C)/SDLauxiliary.h $(S_DIR_C)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_C)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

shadows : $(S_DIR_F)/$(FILE).cpp $(S_DIR_F)/SDLauxiliary.h $(S_DIR_F)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_F)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

dof : $(S_DIR_F)/$(FILE).cpp $(S_DIR_F)/SDLauxiliary.h $(S_DIR_F)/TestModel.h
	$(CC) -D POSTPROCESSING $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_F)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

chromatic_aberration : $(S_DIR_F)/$(FILE).cpp $(S_DIR_F)/SDLauxiliary.h $(S_DIR_F)/TestModel.h
	$(CC) -D POSTPROCESSING -D CHROMATICABERRATION -D CAMERAWARPEFFECT $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_F)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

grain_effect : $(S_DIR_F)/$(FILE).cpp $(S_DIR_F)/SDLauxiliary.h $(S_DIR_F)/TestModel.h
	$(CC) -D POSTPROCESSING -D GRAINEFFECT -D CAMERAWARPEFFECT $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR_F)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)
	$(CC) $(LN_OPTS) -o $(EXEC) $(OBJ) $(SDL_LDFLAGS)

clean:
	rm -f $(B_DIR)/* 
