

function [fmax,best_pos,Convergence_curve]=BOA(n,N_iter,dim,fobj,XX,Temp)

% n is the population size
% N_iter represnets total number of iterations


p=0.8;                       % probabibility switch
power_exponent=0.1;
sensory_modality=0.01;

%Initialize the positions of search agents
Sol=XX;

Fitness=fobj


% Find the current best_pos
[fmax,I]=max(Fitness);
best_pos=Sol(I,:);
S=Sol; 


% Start the iterations -- Butterfly Optimization Algorithm 
for t=1:N_iter,
  display('Interation ',num2str(t));
        for i=1:n, % Loop over all butterflies/solutions
         
          %Calculate fragrance of each butterfly which is correlated with objective function
          Fnew=Get_Functions_details(S(i,:),Temp); ;
          FP=(sensory_modality*(Fnew^power_exponent));   
    
          %Global or local search
          if rand<p,    
            dis = rand * rand * best_pos - Sol(i,:);        %Eq. (2) in paper
            S(i,:)=Sol(i,:)+dis*FP;
           else
              % Find random butterflies in the neighbourhood
              epsilon=rand;
              JK=randperm(n);
              dis=epsilon*epsilon*Sol(JK(1),:)-Sol(JK(2),:);
              S(i,:)=Sol(i,:)+dis*FP;                         %Eq. (3) in paper
          end
          for lk=1:dim
            if(1/(1+exp(-S(i,lk)))>rand)
                S(i,lk)=1;
            else
                S(i,lk)=0;
            end
          end
                    
            % Evaluate new solutions
            Fnew=Get_Functions_details(S(i,:),Temp);  %Fnew represents new fitness values
            
            % If fitness improves (better solutions found), update then
            if (Fnew>Fitness(i)),
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
            end
           
           % Update the current global best_pos
           if Fnew>fmax,
                best_pos=S(i,:);
                fmax=Fnew
           end
         end
            
         Convergence_curve(t,1)=fmax;
         
         %Update sensory_modality
          sensory_modality=sensory_modality_NEW(sensory_modality, N_iter);
        
end


  
function y=sensory_modality_NEW(x,Ngen)
y=x+(0.025/(x*Ngen));



