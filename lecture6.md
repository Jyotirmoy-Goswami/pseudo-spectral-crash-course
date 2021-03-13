# Lecture 6: 
There was an error in the code we were writing in the last lecture :smiley:. Well, I identified the error. 

Look at **Line: 347** in the Subroutine *AdamsBashforth*.
```
    omegak_new(i,j) = omegak(Nh,Ny) + ..........
```
should be replaced as
```
    omegak_new(i,j) = omegak(i,j) + ...........
```
Afterwards, I talked about the basic idea of Kelvin-Helmholtz instability, how to capture it numerically (there are some subtleties here). Finally our baby-code performed well to reproduce the analytical growth rate calculated by P G Drazin and thus passed one very critical test! Then I mentioned some papers, talked about few interesting open problems that can be 'easily' attacked using our simple baby code and finally sketched some aspects of parallel computing that we will be covering in the rest of the course.
