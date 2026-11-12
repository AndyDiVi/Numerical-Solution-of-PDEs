clear;
close all;
clc;

%Impostazione della regione
Px = 12; Py = 8;
Qx = 4; Qy = 2;
Rx = 6; Ry = Qy;

%condizione di Dirichlet sui bordi esterni nord e sud
g0= @(x) 10*(Px-abs(x))^2;

%--------------------------------------------------------------------------
%-----------------------SCELTA DEL PROBLEMA--------------------------------
%--------------------------------------------------------------------------

probl = input(['Quale problema si vuole risolvere? ' ...
    'Immettere 1 per il problema di Laplace, 2 per il problema di Poisson' ...
    ' e 3 per equazione del calore: ']);

while (probl~=1) && (probl~=2) && (probl~=3)
    clc;disp(' -Errore sulla scelta della casistica');
    probl = input(['Quale problema si vuole risolvere? ' ...
        'Immettere 1 per il primo esercizio e 2 per il secondo: ']);
end

if probl == 2
    K0 = 1;
    caso = input(['Quale caso si vuole risolvere? ' ...
        '(Immettere 1 per K discontinuo o 2 per K continuo) ']);
    if caso == 1
        K0 = 110;
    end
    while (caso ~= 1) && (caso ~= 2)
        clc;disp(' -Errore sulla scelta della casistica');
        caso = input(['Quale caso si vuole risolvere? ' ...
            '(Immettere 1 per K discontinuo o 2 per K continuo) ']);
    end
elseif probl==3
    K0 = 1;
    caso = input('Quale metodo si vuole utilizzare? Immettere 1 per EE, 2 per EI, 3 per CN:');
    while (caso~=1) && (caso~=2) && (caso~=3)
        clc;
        disp(' -Errore sulla scelta della casistica');
        caso = input('Quale metodo si vuole utilizzare? Immettere 1 per EE, 2 per EI, 3 per CN:');
    end

end

%--------------------------------------------------------------------------
%-----------------------------GRIGLIA--------------------------------------
%--------------------------------------------------------------------------

% Lettura da input del passo h per la divisione in reticolo
h = input('Inserire un passo di campionamento: ');
% Numero di elementi negli array di divisione
Nx = Px*2/h; Ny = Py*2/h;

% Controllo della correttezza del passo h
if h<1
    %controllo che il numero di punti campionati sia multiplo della
    %dimensione del dominio
    while (mod(2*Px,h)~=0) || (mod(2*Py,h)~=0) || (mod(2,h)~=0) || (mod(4,h)~=0)
        disp('Errore sul passo di campionamento')
        h = input('Inserire un passo di campionamento differente : ');
        Nx = Px*2/h; Ny = Py*2/h;
    end
elseif (h>1) && (h~=2)
    disp('Sono possibili solo passi di ampiezza <=1 o pari a 2');
    disp('Si considera un passo di ampiezza 2')
    h = 2;
    Nx = Px; Ny = Py;
end

if probl==3
    c=0.1; T=10; k=h/2; dim_t=(T/k)+1; lambda=c^2*k/(h^2); %Parametri equazione del calore
    while h^2<4*c^2*k && caso==1
        disp('Errore: con tale scelta non si ha stabilità del metodo')
        h = input('Inserire un passo di campionamento differente: ');
    end
end

%Vettore Larghezza
dim_x = Nx + 2;
x = linspace(-Px-h,Px,dim_x);x=round(x,1);

%Vettore Altezza
dim_y = Ny + 1;
y = linspace(Py,-Py,dim_y);y=round(y,1);

%--------------------------------------------------------------------------
%-----------------------COSTRUZIONE MATRICI G E I--------------------------
%--------------------------------------------------------------------------

%I = sparse(dim_y,dim_x);G = I;
I = zeros(dim_y,dim_x);G = I;

cont = 0; %variabile contatore

if probl ==1   %Problema 1: equazione di Laplace
    for i = 1 : dim_y
        for j = 1 : dim_x

            %Tag punti dentro le regioni gialle
            if ((Qx<abs(x(j))) && (abs(x(j))<Rx) && (abs(y(i))<Qy))
                % Tag dei punti fittizi
                if abs(x(j+1))==Qx || abs(x(j-1))==Qx || abs(x(j-1))==Rx ||  abs(x(j+1))==Rx || y(i+1)==-Ry || y(i-1)==Ry
                    G(i,j)=6;
                    cont = cont+1;
                    I(i,j) = cont;
                % Tag dei punti esclusi dal dominio
                else
                    G(i,j) = 0;
                end

            %Tag dei punti di bordo esterni
            elseif (abs(x(j)) <= Px) && (abs(y(i)) == Py)
                G(i,j) = 2;% Tag dei punti con Diric.Cond g=g0
            elseif (x(j) == Px) && (abs(y(i)) <= Py)
                G(i,j) = 3;% Tag dei punti con Diric.Cond g=4
            elseif (x(j) == -Px) && (abs(y(i)) < Py)
                G(i,j) = 4;% Tag dei punti con Neumann bordo esterno
                cont = cont+1;
                I(i,j) = cont;
            elseif (x(j) == -Px-h) && (abs(y(i)) == Py)
                G(i,j) = 0; % Tag dei punti negli spigoli dai punti fittizi

            % Tag dei Punti con Neumann bordo interno
            elseif ((abs(x(j)) == Qx || abs(x(j)) == Rx) && abs(y(i)) <= Qy) || (abs(y(i)) == Qy && Qx<=abs(x(j)) && abs(x(j))<=Rx)
                %fissa la x o la y e fa variare l'altra coordinata in maniera di beccare tutto il bordo
                G(i,j) = 5;
                cont = cont+1;
                I(i,j) = cont;

            % Tag dei punti fittizi per il calcolo del bordo esterno ovest
            elseif x(j+1) == -Px && abs(y(i)) < Py
                G(i,j) = 6;
                cont = cont+1;
                I(i,j) = cont;

            % Tag dei punti propriamente interni al dominio
            else
                G(i,j) = 1;
                cont = cont+1;
                I(i,j) = cont;

            end
        end
    end

else %Problema 2 e 3: equazioni di Poisson e del calore

    for i = 1 : dim_y
        for j = 1 : dim_x

            %Tag dei punti dentro le regioni gialle
            if ((Qx<abs(x(j))) && (abs(x(j))<Rx) && (abs(y(i))<Qy))
                G(i,j) = 7;
                cont = cont+1;
                I(i,j) = cont;

            %Tag dei punti di bordo interni
            elseif ((abs(x(j)) == Qx || abs(x(j)) == Rx) && abs(y(i)) <= Qy) || (abs(y(i)) == Qy && Qx<=abs(x(j)) && abs(x(j))<=Rx)
                G(i,j)=5;
                cont = cont+1;
                I(i,j) = cont;

            %Tag dei punti di bordo esterni
            elseif (abs(x(j)) <= Px) && (abs(y(i)) == Py)
                G(i,j) = 2;% Tag dei punti con Diric.Cond g=g0
            elseif (x(j) == Px) && (abs(y(i)) <= Py)
                G(i,j) = 3;% Tag dei punti con Diric.Cond g=4
            elseif (x(j) == -Px) && (abs(y(i)) < Py)
                G(i,j) = 4;% Tag dei punti con Neumann bordo esterno
                cont = cont+1;
                I(i,j) = cont;
            elseif (x(j) == -Px-h) && (abs(y(i)) == Py)
                G(i,j) = 0; %Tag dei punti negli spigoli dai punti fittizi

            % Tag dei punti fittizi per il calcolo del bordo esterno ovest
            elseif x(j) == -Px-h && abs(y(i)) < Py
                G(i,j) = 6;
                cont = cont+1;
                I(i,j) = cont;

            else
                G(i,j) = 1;% Tag dei piunti propriamente interni al dominio
                cont = cont+1;
                I(i,j) = cont;

            end
        end
    end

end
clear Px; clear Py; clear Rx; clear Ry; clear Qx; clear Qy;

%%---------Rappresentazione visiva delle matrici---------------------------
figure(); imagesc(G);colorbar; title('Matrice delle geometrie');
%figure(); imagesc(I);colorbar; title('Matrice degli indici');

%Matrice delle geometrie
if probl==1
    figure()
    spy(G==1,'.g')
    hold on
    spy(G==2,'.k')
    hold on
    spy(G==3,'.b')
    hold on
    spy(G==4,'.r')
    hold on
    spy(G==5,'.m')
    hold on
    spy(G==6,'.c')
    hold off
    legend('G=1','G=2','G=3','G=4','G=5', 'G=6')
    title('Matrice delle geometrie')
else
    figure()
    spy(G==1,'.g')
    hold on
    spy(G==2,'.k')
    hold on
    spy(G==3,'.b')
    hold on
    spy(G==4,'.r')
    hold on
    spy(G==5,'.m')
    hold on
    spy(G==6,'.c')
    hold on
    spy(G==7, 'y')
    hold off
    legend('G=1','G=2','G=3','G=4','G=5', 'G=6', 'G=7')
    title('Matrice delle geometrie')
end

%--------------------------------------------------------------------------
%-----------------------COSTRUZIONE SISTEMA LINEARE------------------------
%--------------------------------------------------------------------------

%Inizializzazione di Ah e b per risolvere il sist.lin.
Ah = sparse(cont,cont); b = sparse(cont,1);
%Ah = zeros(cont,cont); b = zeros(cont,1);

if probl == 1
    %Si evita il controllo dei bordi di Dirichlet nord/sud/est, le condizioni
    %vengono implementate a partire dai punti interni
    for i = 2 : dim_y-1
        for j = 1 : dim_x-1

            %punti interni e di bordo, con Neumann
            %in questi punti si applica lo stencil
            if G(i,j)==1 || G(i,j)==4 || G(i,j)==5
                centro=I(i,j);
                sotto=I(i+1,j);
                sopra=I(i-1,j);
                dx=I(i,j+1);
                sx=I(i,j-1);

                Ah(centro,centro)=4;
                Ah(centro,sx)=-1;

                %Gestione dei bordi: assegna a b i termini noti delle
                %condizioni di Dirichlet
                if sopra==0 && G(i-1,j)==2 %Dirichlet nord
                    b(centro)=b(centro)+g0(x(j));
                else
                    Ah(centro,sopra)=-1;
                end
                if sotto==0 && G(i+1,j)==2 %Dirichlet sud
                    b(centro) =b(centro)+g0(x(j));
                else
                    Ah(centro,sotto)=-1;
                end
                if dx==0 && G(i,j+1)==3 %Dirichlet est
                    b(centro) = b(centro)+4;
                else
                    Ah(centro,dx)=-1;
                end

            %Stencil dei punti fittizi
            elseif G(i,j)==6
                centro=I(i,j);

                %Punti fittizi esterni, a sinistra del bordo esterno ovest
                if G(i,j+1)==4
                    dxdx=I(i,j+2);
                    Ah(centro,centro)=1;
                    Ah(centro,dxdx)=-1;

                %Punti fittizi interni
                elseif G(i,j+1)==5 %sul bordo interno destro
                    dxdx=I(i,j+2);
                    Ah(centro,centro)=1;
                    Ah(centro,dxdx)=-1;
                    b(centro)=b(centro)+100*2*h;

                elseif G(i,j-1)==5 %sul bordo interno sinistro
                    sxsx=I(i,j-2);
                    Ah(centro,centro)=1;
                    Ah(centro,sxsx)=-1;
                    b(centro)=b(centro)+100*2*h;

                elseif G(i-1,j)==5 %sul bordo interno sopra
                    susu=I(i-2,j);
                    Ah(centro,centro)=1;
                    Ah(susu,centro)=-1;
                    b(centro)=b(centro)+100*2*h;

                elseif G(i+1,j)==5 %sul bordo interno sotto
                    giugiu=I(i+2,j);
                    Ah(centro,centro)=1;
                    Ah(giugiu,centro)=-1;
                    b(centro)=b(centro)+100*2*h;

                end
            end
        end
    end

else
    %Problema 2 e 3
    for i = 2 : dim_y-1
        for j = 1 : dim_x-1
            centro=I(i,j);
            if x(j)>0
                K=K0;
            else
                K= 1/K0;
            end

            %Punti nella regione gialla: si usa uno stencil diverso da
            %quello dei punti della regione verde
            if G(i,j)==7

                dx = I(i,j+1);
                sx = I(i,j-1);
                sotto = I(i+1,j);
                sopra = I(i-1,j);
                Ah(centro,sopra)=-K;
                Ah(centro,sotto)=-K;
                Ah(centro,dx)=-K;
                Ah(centro,sx)=-K;
                Ah(centro,centro)= 4*K;
                b(centro)=1000*h^2;
            end

            %Stencil punti della regione verde + bordo esterno ovest (Neumann)
            %in questi punti si applica lo stencil del laplaciano
            if G(i,j)==1 || G(i,j)==4
                sotto=I(i+1,j);
                sopra=I(i-1,j);
                dx=I(i,j+1);
                sx=I(i,j-1);

                Ah(centro,centro)=4;
                Ah(centro,sx)=-1;
                
                %Gestione dei bordi: assegna a b i termini noti delle
                %condizioni di Dirichlet
                if G(i-1,j)==2 %dirichlet nord
                    b(centro)=b(centro)+g0(x(j));
                else
                    Ah(centro,sopra)=-1;
                end

                if G(i+1,j)==2 %dirichlet sud
                    b(centro) =b(centro)+g0(x(j));
                else
                    Ah(centro,sotto)=-1;
                end

                if G(i,j+1)==3 %dirichlet est
                    b(centro) = b(centro)+4;
                else
                    Ah(centro,dx)=-1;
                end
            end

            %Regione gialla
            %punti di bordo interni
            if G(i,j)==5
                dx = I(i,j+1);
                sx = I(i,j-1);
                sotto = I(i+1,j);
                sopra = I(i-1,j);

                if G(i,j-1)==1
                    Ah(centro,dx)=-K;
                    Ah(centro,sx)=-1;

                elseif G(i,j+1)==1
                    Ah(centro,sx)=-K;
                    Ah(centro,dx)=-1;

                elseif G(i-1,j)==1
                    Ah(centro,sotto)=-K;
                    Ah(centro,sopra)=-1;

                elseif G(i+1,j)==1
                    Ah(centro,sopra)=-K;
                    Ah(centro,sotto)=-1;
                end

                Ah(centro,centro)= K+1; 
            end

            %punti fittizi
            %non va applicato lo stencil del laplaciano
            if G(i,j)==6
                centro=I(i,j);
                dxdx=I(i,j+2);
                Ah(centro,centro)=1;
                Ah(centro,dxdx)=-1;
            end
        end
    end

end
clear K; clear K0;
clear sopra; clear sotto; clear centro; clear sx; clear dx;
clear sxsx; clear dxdx; clear susu; clear giugiu;
% %------------Verifica Correttezza------------------------------------------
%if abs(det(Ah))>0
%     disp('Tutto corretto');
% else
%     disp('e niente, avete sbagliato qualcosa...');
%     disp('La matrice Ah ha rango '); disp(rank(Ah));
%     z = sum(abs(Ah'));
%     ind = find(z==0);
%     if isempty(ind)
%         disp('e non ha righe nulle!');
%     else
%         disp('e le righe incriminate sono contenute nella variabile "ind" ');
%     end
% end

%%---------------Rappresentazione visiva zoom di Ah------------------------
% figure; imagesc(Ah);colorbar;title('Matrice Ah');
%
% %Plot zoom matrice Ah
% if probl==1
%     A=Ah(500:650, 500:650);
%     figure()
%     spy(A==-1,'.k')
%     hold on
%     spy(A==1,'.r')
%     hold on
%     spy(A==4,'.b')
%     hold off
%     legend({'Ah =-1', 'Ah = 1', 'Ah = 4'}, 'Interpreter', 'latex')
% else
%     if caso==2
%         A=Ah(500:650, 500:650);
%         figure()
%         spy(A==2,'.g')
%         hold on
%         spy(A==-1,'.k')
%         hold on
%         spy(A==1,'.r')
%         hold on
%         spy(A==4,'.b')
%         hold off
%         legend({'Ah=2', 'Ah = -1', 'Ah = 1', 'Ah = 4'}, 'Interpreter', 'latex')
%     else
%         A=Ah(500:650, 500:650);
%         figure()
%         spy(A~=0,'.g')
%         hold on
%         spy(A==-1,'.k')
%         hold on
%         spy(A==1,'.r')
%         hold on
%         spy(A==4,'.b')
%         hold off
%         legend({'Ah $\neq$ 0', 'Ah = -1', 'Ah = 1', 'Ah = 4'}, 'Interpreter', 'latex')
%     end
% end
% xticks(1:10:151) 
% xticklabels(500:10:650) 
% yticks(1:10:151) 
% yticklabels(500:10:650)

%-------------------------Analisi matrice Ah-------------------------------
if probl==1 || probl==2
    %Singolarità ed elementi distinti
    Ah_elements=unique(full(Ah));
    fprintf('La matrice Ah è non singolare con determinante %s .\n',det(full(Ah)));
    disp('Gli elementi distinti sono: ');disp(Ah_elements); clear Ah_elements;

    %Sparsità Ah
    non_null = 0;
    for i=1:size(Ah,1)
        for j=1:size(Ah,2)
            if(Ah(i,j) ~= 0)
                non_null = non_null + 1;
            end
        end
    end
    aux=(cont^2-non_null)/cont^2; %sparsità
    %Si nomina "aux" l'output poichè utile solo in questo caso, verrà poi
    %utilizzato lo stesso nome in seguito per output poco utilizzati,
    %sovrascrivendo la variabile e allocando meno memoria
    fprintf('La matrice Ah è sparsa al %s percento .\n',aux*100); clear non_null;
end
clear cont;
%--------------------------------------------------------------------------
%------------------------------SOLUZIONE-----------------------------------
%--------------------------------------------------------------------------
%Risoluzione sistema lineare: metodo diretto
tic
u0=Ah\b;
toc

%Verifica condizioni: metodi iterativi
%%raggio spettrale: (da utilizzare solo se Ah inizializzata come zeros)
%rho=max(eig(Ah));
%%simmetria della matrice:
%isequal(Ah, Ah')

%Soluzione = sparse(size(I));
Soluzione = zeros(size(I));
for i = 1:dim_y
    for j = 1:dim_x
        if I(i,j)>0 && G(i,j)~=6
            Soluzione(i,j) = u0(I(i,j));
        elseif G(i,j) == 2
            Soluzione(i,j) = g0(x(j));

        elseif  G(i,j)==3
            Soluzione(i,j) = 4;
        else
            Soluzione(i,j) =NaN;

        end
    end
end

%--------------------------------------------------------------------------
%--------------------------CONFRONTO SOLUZIONI FEM/FDM---------------------
%----------------------------ESERCIZIO 3-----------------------------------
if probl==1 || probl==2
%     %Importazione dei risultati FEM: soluzioni u e geometrie g al variare di h
%     if probl==1
%         FEM_mesh=importdata('FEM solutions\g1.mat');
%         if h==1
%             u=importdata('FEM solutions\u1_h1.mat');
%         elseif h==0.5
%             u=importdata('FEM solutions\u1_h05.mat');
%         elseif h==0.2
%             u=importdata('FEM solutions\u1_h02.mat');
%         end
%     else
%         if probl==2
%             if caso==1
%                 FEM_mesh=importdata('FEM solutions\g2a.mat');
%                 if h==1
%                     u=importdata('FEM solutions\u2a_h1.mat');
%                 elseif h==0.5
%                     u=importdata('FEM solutions\u2a_h05.mat');
%                 elseif h==0.2
%                     u=importdata('FEM solutions\u2a_h02.mat');
%                 end
% 
%             else
%                 FEM_mesh=importdata('FEM solutions\g2b.mat');
%                 if h==1
%                     u=importdata('FEM solutions\u2b_h1.mat');
%                 elseif h==0.5
%                     u=importdata('FEM solutions\u2b_h05.mat');
%                 elseif h==0.2
%                     u=importdata('FEM solutions\u2b_h02.mat');
%                 end
%             end
%         end
%     end
% 
%     %Importazione dei dati della triangolazione
%     [p,~,t]=initmesh(FEM_mesh, 'Hmax', h);
% 
%     %Selezione di punti test su cui effettuare il confronto delle soluzioni
%     punti = [0,-12,-6,-4,4,6;
%              0, 0,  0, 0,0,0];
% 
%     F = pdeInterpolant(p,t,u);
%     ris_FEM= evaluate(F,punti);
%     clear p; clear t; clear punti; clear F;
% 
%     C_0_x = ((Nx+2)/2)+1; %coordinata x dello zero
%     C_0_y = (Ny/2)+1; %coordinata y dello zero
% 
%     %Vettore delle soluzioni FDM nei punti
%     ris_FDM = [Soluzione(C_0_y,C_0_x);Soluzione(C_0_y,C_0_x - 12/h);
%                Soluzione(C_0_y,C_0_x -6/h);Soluzione(C_0_y,C_0_x-4/h);
%                Soluzione(C_0_y,C_0_x+4/h);Soluzione(C_0_y,C_0_x+6/h)];
%     clear C_0_x; clear C_0_y;
% 
%     %Errore assoluto: differenza tra le due soluzioni (FEM e FDM)
%     Diff_Sol = abs(ris_FEM-ris_FDM);
%     %Errore relativo percentuale
%     aux= (ris_FEM+ris_FDM)/2; %risultato medio, denominato aux per allocare meno memoria
%     DIff_Sol_Perc = Diff_Sol./aux*100;
    %--------------------------------------------------------------------------
    %--------------------------ANALISI RISULTATI-------------------------------
    %--------------------------------------------------------------------------

    figure;imagesc(Soluzione);colorbar;title(['Soluzione imagesc con h=', num2str(h)], 'FontSize', 35);
    figure;contourf(Soluzione);colorbar;title(['Soluzione contour con h=', num2str(h)], 'FontSize', 35);
    
    %Range soluzione: massimo e punti di massimo
    aux=max(Soluzione, [], 'all'); %valore massimo, denominato aux per allocare meno memoria
    [max_y, max_x] = find(Soluzione==aux);

    fprintf('La soluzione assume valore massimo %d nei punti di coordinate:\n', aux);
    for i = 1:size(max_x)
        fprintf('x: %d, y: %d\n', x(max_x(i)), y(max_y(i)));
    end
    clear max_y; clear max_x;
    
    %Range soluzione: massimo e punti di massimo
    aux=min(Soluzione, [], 'all'); %valore minimo
    [min_y, min_x] = find(Soluzione==aux);
    
    fprintf('La soluzione assume valore minimo %d nei punti di coordinate:\n', aux);
    for i = 1:size(min_x)
        fprintf('x: %d, y: %d\n', x(min_x(i)), y(min_y(i)));
    end
    clear min_y; clear min_x;

    % Plot soluzione
    figure
    [X,Y] = meshgrid(x,y);
    mesh(X,Y,Soluzione);
    colormap('parula');colorbar;
    xlabel('x');ylabel('y');title(['Soluzione 3D con h=', num2str(h)], 'FontSize', 35 );
    view(-22,33)

    %Plot antigradiente
    figure
    [fx,fy] = gradient(Soluzione);
    contour(X,Y,Soluzione,'ShowText','on');c.LineWidth = 2;
    hold on
    quiver(X,Y,-fx,-fy) %E=-gradiente u
    hold off
    colormap('parula');colorbar;
    xlabel('x');ylabel('y');title(['Anti-gradiente con h=', num2str(h)], 'FontSize', 35);

    %Modulo del gradiente
    mod_grad =zeros(size(G,1), size(G,2));
    for i=1:size(G,1)
        for j=1:size(G,2)
            mod_grad(i,j) = sqrt(fx(i,j)^2+fy(i,j)^2);
        end
    end
    clear fx; clear fy;

    %Massimo del modulo del gradiente
    aux= max(mod_grad, [], 'all');
    [y_grad, x_grad] = find(mod_grad==aux);
    fprintf('Il gradiente assume valore massimo in modulo %d nei punti di coordinate:\n', aux);
    for i = 1:size(x_grad)
        fprintf('x: %d, y: %d\n', x(x_grad(i)), y(y_grad(i)));
    end

    %Minimo del modulo del gradiente
    aux= min(mod_grad, [], 'all');
    [y_grad, x_grad] = find(mod_grad==aux);
    fprintf('Il gradiente assume valore minimo in modulo %d nei punti di coordinate:\n', aux);
    for i = 1:size(x_grad)
        fprintf('x: %d, y: %d\n', x(x_grad(i)), y(y_grad(i)));
    end
    clear y_grad; clear x_grad;clear aux;

%     %Plot del modulo del gradiente
%     figure
%     [X,Y] = meshgrid(x,y);
%     mesh(X,Y,mod_grad)
%     colormap('parula');colorbar;
%     xlabel('x');ylabel('y');title(['Modulo del gradiente con h=', num2str(h)], 'FontSize', 35);
%     view(15,65)
end
%--------------------------------------------------------------------------
%--------------------------------ESERCIZIO 4-------------------------------
%----------------------------EQUAZIONE DEL CALORE--------------------------
%-------------------------------------------------------------------------
if probl==3
    A = zeros(length(b)); b = zeros(size(b)); %stencil per il tempo n+1

    Sol=zeros(length(b),dim_t);
    Sol(:,1)=u0;

    if caso==1  %Eulero esplicito
        for i = 2 : dim_y-1
            for j = 1 : dim_x-1

                centro=I(i,j);
                A(centro,centro)=1-(4*lambda);

                %punti interni bordo esterno ovest (Neumann)
                if G(i,j)==1 || G(i,j)==4 || G(i,j)==5 || G(i,j)==7
                    sotto=I(i+1,j);
                    sopra=I(i-1,j);
                    dx=I(i,j+1);
                    sx=I(i,j-1);


                    A(centro,sx) = lambda;

                    %Gestione dei bordi: assegna a b i termini noti delle
                    %condizioni di Dirichlet
                    if G(i-1,j)==2 %dirichlet nord
                        b(centro)=g0(x(j))*lambda;

                    else
                        A(centro,sopra) = lambda;

                    end

                    if G(i+1,j)==2 %dirichlet sud
                        b(centro)=g0(x(j))*lambda;
                    else
                        A(centro,sotto) = lambda;
                    end

                    if G(i,j+1)==3 %dirichlet est
                        b(centro)=4*lambda;
                    else
                        A(centro,dx) = lambda;
                    end
                end
            end
        end
        %Sistema lineare Eulero esplicito
        for j=1:dim_t-1
            Sol(:,j+1)=(A*Sol(:,j))+b;
        end

    elseif caso==2 %Eulero Implicito
        for i = 2 : dim_y-1
            for j = 1 : dim_x-1

                centro=I(i,j);
                A(centro,centro)=1+(4*lambda);

                %punti interni bordo esterno ovest (Neumann)
                if G(i,j)==1 || G(i,j)==4 || G(i,j)==5 || G(i,j)==7
                    sotto=I(i+1,j);
                    sopra=I(i-1,j);
                    dx=I(i,j+1);
                    sx=I(i,j-1);


                    A(centro,sx) = -lambda;

                    %Gestione dei bordi: assegna a b i termini noti delle
                    %condizioni di Dirichlet
                    if G(i-1,j)==2 %dirichlet nord
                        b(centro)=g0(x(j))*lambda;

                    else
                        A(centro,sopra) = -lambda;

                    end

                    if G(i+1,j)==2 %dirichlet sud
                        b(centro)=g0(x(j))*lambda;
                    else
                        A(centro,sotto) = -lambda;
                    end

                    if G(i,j+1)==3 %dirichlet est
                        b(centro)=4*lambda;
                    else
                        A(centro,dx) = -lambda;
                    end
                end
            end
        end
        
        %Sistema lineare Eulero implicito
        for j=1:dim_t-1
            Sol(:,j+1)=A\(Sol(:,j)+b);
        end
    else %Crank Nicolson
        B = A; bb = zeros(size(b)); %stencil per il tempo n
        for i = 2 : dim_y-1
            for j = 1 : dim_x-1

                centro=I(i,j);
                A(centro,centro)=1+(2*lambda);
                B(centro,centro)=1-(2*lambda);

                %punti interni bordo esterno ovest (Neumann)
                if G(i,j)==1 || G(i,j)==4 || G(i,j)==5 || G(i,j)==7
                    sotto=I(i+1,j);
                    sopra=I(i-1,j);
                    dx=I(i,j+1);
                    sx=I(i,j-1);


                    A(centro,sx)=-(lambda/2);
                    B(centro,sx)=(lambda/2);

                    %Gestione dei bordi: assegna a b i termini noti delle
                    %condizioni di Dirichlet
                    if G(i-1,j)==2 %dirichlet nord
                        b(centro)=-(g0(x(j))*lambda/2);
                        bb(centro)=(g0(x(j))*lambda/2); 
                    else
                        A(centro,sopra)=-(lambda/2);
                        B(centro,sopra)=(lambda/2);
                    end

                    if G(i+1,j)==2 %dirichlet sud
                        b(centro)=-(g0(x(j))*lambda/2);
                        bb(centro)=(g0(x(j))*lambda/2);
                    else
                        A(centro,sotto)=-(lambda/2);
                        B(centro,sotto)=(lambda/2);
                    end

                    if G(i,j+1)==3 %dirichlet est
                        b(centro)=-(2*lambda);
                        bb(centro)=+(2*lambda);
                    else
                        A(centro,dx)=-(lambda/2);
                        B(centro,dx)=(lambda/2);
                    end
                end
            end
        end

        %Sistema lineare Crank-Nicolson
        for j=1:dim_t-1
            termine_noto = B*Sol(:,j)+bb;
            termine_noto = termine_noto - b;
            Sol(:,j+1)=A\termine_noto;
        end
    end
    clear centro; clear dx; clear sx; clear sopra; clear sotto;
%----------------------SOLUZIONE------------------------------------------
    Soluzione = zeros(size(I,1),size(I,2),dim_t);
    for z=1:dim_t
        for i = 1:dim_y
            for j = 1:dim_x
                if (I(i,j)>0) && (G(i,j)~=6)
                    Soluzione(i,j,z) = Sol(I(i,j),z);
                elseif G(i,j) == 2
                    Soluzione(i,j,z) = g0(x(j));

                elseif G(i,j) == 3
                    Soluzione(i,j,z) = 4;

                else
                    Soluzione(i,j,z) = NaN;
                end
            end
        end
    end
    clear g0;

%-----------------------GRAFICI/PLOT---------------------------------------
    %Calcoliamo il max per aggiustare le scale dei plot
    max=max(Soluzione(:,:,1), [], 'all'); 
    [X,Y] = meshgrid(x,y);

    %N.B. vogliamo il codice funzionante per ogni eventuale tempo T così da
    %testare diversi tempi, per questo motivo non ci limitiamo a mettere
    %valori fissati

    %Plot a tempo successivo a quello iniziale
    figure(); %2D
    imagesc(Soluzione(:,:,round(2/k)+1), [0 max]);colorbar;title(sprintf('Soluzione a tempo t=%d', 1), 'FontSize', 35);

    figure; %3D
    mesh(X,Y,Soluzione(:,:,round(2/k)+1))
    colormap('parula');
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Soluzione a tempo t=%d', 1), 'FontSize', 35);
    clim([0 max]);

    %Plot a metà tempo
    figure(); %2D
    imagesc(Soluzione(:,:,round(dim_t/2)+1), [0 max]);colorbar;title(sprintf('Soluzione a tempo t=%d', round(k*(dim_t/2-1))), 'FontSize', 35);
    figure; %3D
    mesh(X,Y,Soluzione(:,:,round(dim_t/2)+1))
    colormap('parula');
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Soluzione a tempo t=%d', round(k*(dim_t/2-1)) ),'FontSize', 35);
    clim([0 max]);

   %Plot a tempo finale
    figure(); %2D
    imagesc(Soluzione(:,:,dim_t),[0 max]);colorbar;title(sprintf('Soluzione a tempo t=%d', T), 'FontSize', 35);

    figure; %3D
    mesh(X,Y,Soluzione(:,:,dim_t))
    colormap('parula');
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Soluzione a tempo t=%d', T), 'FontSize', 35);
    clim([0 max]);
    
%%-------------------------------------------------------------------------
%%-----------------------------------VIDEO 2D-------------------------------
%%-------------------------------------------------------------------------
%         videoFileName = 'EquazioneCalore.mp4';
%         video = VideoWriter(videoFileName,'MPEG-4');
%         frameRate = 10;
%         video.FrameRate = frameRate;
%         open(video);
%         for z = 1:dim_t
%             figure();
%             imagesc(Soluzione(:,:,z), [0 max]);
%              colorbar;
%             title(sprintf('Soluzione a tempo t=%.2f', round((z-1)*k,3)));
%             frame = getframe(gcf);
%             writeVideo(video, frame);
%             close;
%         end
%         close(video);
%         clear videoFileName; clear video; clear frame; clear frameRate;
%     
%     
%%--------------------------------------------------------------------------
%%----------------------------------VIDEO 3D--------------------------------
%%-------------------------------------------------------------------------
%         videoFileName = 'EquazioneCalore3D.mp4';
%         video = VideoWriter(videoFileName,'MPEG-4');
%         frameRate = 10;
%         video.FrameRate = frameRate;
%         open(video);
%     
%         for z = 1:dim_t
%             figure();
%             mesh(X,Y,Soluzione(:,:,z))
%             colormap('parula');
%             colorbar;
%             xlabel('x');
%             ylabel('y');
%             title(sprintf('Soluzione a tempo t=%.2f', round(k*(z-1),3)));
%             clim([0 max]);
%             frame = getframe(gcf);
%             writeVideo(video, frame);
%             close;
%         end
%         close(video);
%         clear videoFileName; clear video; clear frame; clear frameRate;
%%%-------------------------------------------------------------------------
%%-------------------------------MAX E MIN PARABOLICA----------------------
%%-------------------------------------------------------------------------
%         %ANALISI MASSIMI

%         %Tempo iniziale
%         aux=max(max(Soluzione(:,:,1)));
%         %aux=max(Soluzione(:,:,1), [], 'all'); %valore massimo, denominato aux per allocare meno memoria
%         [max_y, max_x] = find(Soluzione(:,:,1)==aux);
%         
%         fprintf('La soluzione assume valore massimo %d nei punti di coordinate:\n', aux);
%         for i = 1:size(max_x)
%             fprintf('x: %d, y: %d\n', x(max_x(i)), y(max_y(i)));
%         end

%         %Tempo 2
%         aux=max(max(Soluzione(:,:,round(2/k)+1)));
%         %aux=max(Soluzione(:,:,(2/k)+1), [], 'all'); %valore massimo, denominato aux per allocare meno memoria
%         [max_y, max_x] = find(Soluzione(:,:,(2/k)+1)==aux);
%         
%         fprintf('La soluzione assume valore massimo %d nei punti di coordinate:\n', aux);
%         for i = 1:size(max_x)
%             fprintf('x: %d, y: %d\n', x(max_x(i)), y(max_y(i)));
%         end

%         %Tempo 5
%         aux=max(max(Soluzione(:,:,round(5/k)+1)));
%         %aux=max(Soluzione(:,:,(5/k)+1), [], 'all'); %valore massimo, denominato aux per allocare meno memoria
%         [max_y, max_x] = find(Soluzione(:,:,(5/k)+1)==aux);
%         
%         fprintf('La soluzione assume valore massimo %d nei punti di coordinate:\n', aux);
%         for i = 1:size(max_x)
%             fprintf('x: %d, y: %d\n', x(max_x(i)), y(max_y(i)));
%         end

%         %Tempo 10
%         aux=max(max(Soluzione(:,:,dim_t)));
%         %aux=max(Soluzione(:,:,dim_t), [], 'all'); %valore massimo, denominato aux per allocare meno memoria
%         [max_y, max_x] = find(Soluzione(:,:,dim_t)==aux);
%         
%         fprintf('La soluzione assume valore massimo %d nei punti di coordinate:\n', aux);
%         for i = 1:size(max_x)
%             fprintf('x: %d, y: %d\n', x(max_x(i)), y(max_y(i)));
%         end
% 
%         %ANALISI MINIMI

%         %Tempo iniziale
%         aux=min(min(Soluzione(:,:,1)));
%         %aux=min(Soluzione(:,:,1), [], 'all'); %valore minimo, denominato aux per allocare meno memoria
%         [min_y, min_x] = find(Soluzione(:,:,1)==aux);
%         
%         fprintf('La soluzione assume valore min %d nei punti di coordinate:\n', aux);
%         for i = 1:size(min_x)
%             fprintf('x: %d, y: %d\n', x(min_x(i)), y(min_y(i)));
%         end

%         %Tempo 2
%         aux=min(min(Soluzione(:,:,round(2/k)+1)));
%         %aux=min(Soluzione(:,:,(2/k)+1), [], 'all'); %valore minimo, denominato aux per allocare meno memoria
%         [min_y, min_x] = find(Soluzione(:,:,(2/k)+1)==aux);
%         
%         fprintf('La soluzione assume valore min %d nei punti di coordinate:\n', aux);
%         for i = 1:size(min_x)
%             fprintf('x: %d, y: %d\n', x(min_x(i)), y(min_y(i)));
%         end

%         %Tempo 5
%         aux=min(min(Soluzione(:,:,round(5/k)+1)));
%         %aux=min(Soluzione(:,:,(5/k)+1), [], 'all'); %valore minimo, denominato aux per allocare meno memoria
%         [min_y, min_x] = find(Soluzione(:,:,(5/k)+1)==aux);
%         
%         fprintf('La soluzione assume valore min %d nei punti di coordinate:\n', aux);
%         for i = 1:size(min_x)
%             fprintf('x: %d, y: %d\n', x(min_x(i)), y(min_y(i)));
%         end

%         %Tempo 10
%         aux=min(min(Soluzione(:,:,dim_t)));
%         %aux=min(Soluzione(:,:,dim_t), [], 'all'); %valore minimo, denominato aux per allocare meno memoria
%         [min_y, min_x] = find(Soluzione(:,:,dim_t)==aux);
%         
%         fprintf('La soluzione assume valore min %d nei punti di coordinate:\n', aux);
%         for i = 1:size(min_x)
%             fprintf('x: %d, y: %d\n', x(min_x(i)), y(min_y(i)));
%         end

end
clear i; clear j; clear caso; clear probl; clear max;