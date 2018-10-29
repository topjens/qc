%{
#include <stdio.h>
#include <stdlib.h>
%}
                        
%union
{
    int i;
    double f;
    char *s;
}
                        
%token<i> INTEGER SHELL
%token<f> FLOAT
%token<s> ELEMENT

%start file

%%
file:           /* nothing */
        |       elements
        ;

elements:       element
        |       element elements
        ;

element:        ELEMENT shells
                {
                    printf("New element %s!\n", $1);
                }
        ;

shells:         shell
        |       shell shells
        ;

shell:          SHELL INTEGER shell_details
        ;

shell_details:  shell_detail
        |       shell_detail shell_details
        ;

shell_detail:   INTEGER FLOAT FLOAT
        ;

%%
int main(int argc, char *argv[])
{
	yyparse();
}


