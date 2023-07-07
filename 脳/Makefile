MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
CC = g++
CFLAGS = -O2 -Wall -std=c++17 -Wno-unknown-pragmas

.PHONY: 脳 clean


%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

脳: json/jsoncpp.o common.o 脳.o
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -rf *.o
