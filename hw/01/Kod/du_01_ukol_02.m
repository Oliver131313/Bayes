%% __________________________________________Domaca uloha, ukol 2 Overview______________________________________________
% 
% (1)   Nastavenie potrebnych premennych 
%
% (2)   Vytvorenie vsetkych kombinacii golov s obmedzenim na maximalny pocet golov
%   (i)     Vsetky mozne kombinacie 
%   (ii)    Obmedzenie iba na kombinacie, ktorych sucet golov je <= hodnote stanovenej v max_goly
% 
% (3)   Samotna simulacia
%   (i)     Vygenerovanie vysledku (1=Barca; 2=Remiza; 3=Real) s
%           respektovanim danych pravdepodobnosti (p(1) = 0.4; p(2) = 0.25; p(3) = 0.35) 
%           + Vygenerovanie intervalu celkoveho poctu golov (<3; >=3) s respektovanim pravdepodobnosti 
%             (p(<3) = 0.5; p(>=3) = 0.5)
%   (ii)    Simulovacie podla toho aky nastal vysledok
%       (a)     Vysledok = 1 (Barca) 
%                   --> Vypocet ObaGolMnozina a SkoreMnozina pre pocet golov <3 a >=3
%       (b)     Vysledok = 2 (Remiza)
%                   --> Vypocet ObaGolMnozina a SkoreMnozina pre pocet golov <3 a >=3
%       (c)     Vysledok = 3 (Real)
%                   --> Vypocet ObaGolMnozina a SkoreMnozina pre pocet golov <3 a >=3
% 
% (4) Vysledky a vizualizacia 
clear ALL;
clear variables;
%% 1. Nastavenie premennych
p_barca = 0.4;
p_remiza = 0.25;
p_real = 0.35;
p_menej_ako_3 = 0.5;
p_3_a_viac = 0.5;
S = 100000;         % Pocet simulacii
max_goly = 10;      % Stanovenie maximalneho poctu golov

%% 2. Vygenerovanie vsetkych moznych kombinacii vysledkov
v1 = [];
v2 = [];
for a = 0:max_goly
    for b = 0:max_goly
        v1(:, end + 1) = a;
        v2(:, end + 1) = b;
    end
end

% Vsetky mozne vysledky ak by mohlo padat aj VIAC ako max_goly
V = [v1', v2'];
% Pridanie stlpaca pre sucet golov
V = [V, sum(V, 2)];
% Filter takych zapasov kde je sucet golov <= 10
filter = V(:, 3) <= max_goly;
combs = V(filter, :);

%% 3. For loop - simulacia
P = 0;                      
r = 0;
% Samotna simulacia
for s = 1:S
    % Vygeneruj vysledok zapasu
    r_zapas = rand;
    prob_zapas = [p_barca, p_remiza, p_real];                       % Vektor pravdepodobnosti vysledkov
	vysledok = sum(r_zapas >= cumsum([0, prob_zapas]));             % 1 = vyhra Barca; 2 = Remizu; 3 = Real
    % Vygeneruj ci bol sucet golov >= 3 alebo < 3
    r_goly = rand;
    prob_goly = [p_menej_ako_3, p_3_a_viac];                        % Vektor pravdepodobnosti suctu golov
    poc_golov = sum(r_goly >= cumsum([0, prob_goly]));
    
    % Ak vyhra barcelona
    if vysledok == 1
        % Ak je pocet golov vacsi alebo rovny 3
        if poc_golov == 2
            r = r+1;
            % Filter vysledkov kedy vyhra Barca a sucet golov >= 3
            filter_1 = (combs(:, 1) > combs(:, 2)) & (combs(:, 3) >= 3);  
            combs_1 = combs(filter_1, 2);
            p_Rgoly = sum(combs_1 > 0) / size(combs_1, 1);          % Proporcia zapasov kedy dal Real aspon 1 gol
            P = P + p_Rgoly;                                        % Pravdepodobnost pozadovaneho stavu

        end
        
    % Ak bude remiza        
    elseif vysledok == 2
        % Ak je pocet golov < 3
        if poc_golov == 1
            r = r+1;
            % Filter: Remizy a zaroven sucet golov < 3
            filter_2 = (combs(:, 1) == combs(:, 2)) & (combs(:, 3) < 3);
            combs_2 = combs(filter_2, 3);
            p_RemizaGoly = sum(combs_2 == 2) / size(combs_2, 1);     % Proporcia zapasov, kedy oba timy skorovali  
            P = P + p_RemizaGoly;                                    % Pravdepodobnost pozadovaneho stavu
        % Ak je pocet golov >= 3 a nastala remiza, tak je jasne, ze oba
        % timy museli skorovat
        elseif poc_golov == 2
            r = r+1;
            % Filter: Remizy a zaroven je sucet golov >= 3
            filter_2 = (combs(:, 1) == combs(:, 2)) & (combs(:, 3) >= 3);
            combs_2 = combs(filter_2, 3);
            p_RemizaGoly = sum(combs_2 >= 3) / size(combs_2, 1);     % Proporcia zapasov, kedy oba timy skorovali (= 1)
            P = P + p_RemizaGoly;                                    % Pravdepodobnost pozadovaneho stavu
        end
    % Ak vyhra Real
    elseif vysledok == 3
        % Ak je pocet golov < 3
        if poc_golov == 2
            r = r+1;
            ObaGolMnozina = ObaGolMnozina + 1;
            % Filter: Vyhra Real Madrid a pocet golov >= 3
            filter_3 = (combs(:, 1) < combs(:, 2)) & (combs(:, 3) >= 3);
            combs_3 = combs(filter_3, 1);
            % Minimalny a maximalny pocet golov, ktore mohla Barca vstrelit
            % za splnenia predoslych podmienok
            p_Bgoly = sum(combs_3 > 0) / size(combs_3, 1);            % Proporcia zapasov, kedy oba timy skorovali
            P = P + p_Bgoly;                                          % Pravdepodobnost pozadovaneho stavu
         end
    else 
        disp("Niekde mas chybu, pretoze by si mal vycerpat vsetky mozne stavy")
    end
end

%% 4. Vysledky
E_ObaGoly = P / r                                            % Spravodliva pravdepodobnost (ocakavana hodnota)
Kurz_ObaGoly = 1/E_ObaGoly;                                  % Spravodlivy kurz stavkovej kancelarie
E_NieObaGoly = 1 - E_ObaGoly;
Kurz_NieObaGoly = 1/E_NieObaGoly;

fprintf('\nPravdepodobnost, ze oba timy v zapase skoruju:\n\t%.3f', ... 
        E_ObaGoly)
fprintf('\nSpravodlivy kurz, ktory by na tuto udalost mala vypisat stavkova kancelaria:\n\t%.2f\n',... 
        Kurz_ObaGoly)
fprintf('\nPravdepodobnost, ze oba timy v zapasae NEskoruju:\n\t%.3f', ... 
        E_NieObaGoly)
fprintf('\nSpravodlivy kurz, ktory by na tuto udalost mala vypisat stavkova kancelaria:\n\t%.2f\n',... 
        Kurz_NieObaGoly)