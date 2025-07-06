#!/bin/bash
# The above line is called a shebang and it tells the system to use bash to interpret this script   


# Printing:
# echo is like print in Python
echo -e "Hello, World\n"

# Prompting:
# read lets you take input from the user
# -p allows you to prompt the user for input
# read -p "Please enter your name  : " name
# echo ""
# read -p "Please enter your age  : " age
# echo ""
# read -p "Please enter your sex. Male/Female  : " sex
# echo ""
# echo "So you're a $age year old $sex called $name"

# if Statements:
echo "Please enter type of fruit"
read fruit
# $ indicates a variable
if [ $fruit = apple ]
        then echo "Good, I like Apples"
elif [ $fruit = pear ]
        then echo "Good, I like Pears"
elif [ $fruit = banana ]
        then echo "Good, I like Bananas"
        else echo "Oh no, I hate Oranges!"
fi
# fi = finish

# Creating functions

