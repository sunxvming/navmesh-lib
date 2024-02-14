.SUFFIXES: .cpp .o

CC=$(CROSS)clang $(TARGET_FLAGS)  -lstdc++
AR=$(CROSS)ar -cr
SRCS=Circle.cpp Edge.cpp Math.cpp NavPath.cpp Path.cpp Point.cpp Polygon.cpp Triangle.cpp 


OBJS=$(SRCS:.cpp=.o)
EXEC=libNavPath.a


start: $(OBJS)
	$(AR)  $(EXEC) $(OBJS)
	@echo '----------------ok------------$(CC)'

.c.o:
	$(CC) -Wall -g -o $@ -c $<
clean:
	rm -f $(OBJS)
	rm -f libNavPath.a
