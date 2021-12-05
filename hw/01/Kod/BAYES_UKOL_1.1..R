s = 100000         # pocet simulaci
r = 0              # inicializace pocitadla zadaneho stavu 
E_Bwin <- 0        # inicializace pocitadla pravdepodobnosti vyhry boba ... pozdeji zprumerovana 

CIL_HRY <- as.numeric(readline(prompt= "Zadejte pocet bodu nutny k vitezstvi:"))  
BODY_ALICE <- as.numeric(readline(prompt= "Zadejte pocet bodu Alice:"))
BODY_BOB <- as.numeric(readline(prompt= "Zadejte pocet bodu Boba:"))
MIN_POCET_HER_WIN_BOB <- CIL_HRY - BODY_BOB                           # minimalni pocet her potrebny k vyhre boba od zastaveni hry
MAX_POCET_HER_WIN_BOB <- (CIL_HRY-1)*2 - BODY_BOB - BODY_ALICE + 1    # maximalni pocet her potrebny k vyhre boba od zastaveni hry

for (i in 1:s){
  p <- runif(1, min = 0, max = 1)
  y <- runif(BODY_ALICE + BODY_BOB, min = 0, max = 1)                 # BODY_ALICE + BODY_BOB = pocet her odehranych v moment zastaveni hry 
  if (sum(y < p) == BODY_ALICE) {
    r = r+1
    for (k in MIN_POCET_HER_WIN_BOB:MAX_POCET_HER_WIN_BOB){           # cyklus pricitajici pravdepodobnosti vyhry boba za zadaneho stavu - cyklus probehne pro kazdy pocet her od zastaveni hry, ktery muze vest k vyhre boba
      E_Bwin = E_Bwin + (1-p)**(MIN_POCET_HER_WIN_BOB) * (p)**(k-MIN_POCET_HER_WIN_BOB) * choose(ifelse(k == min(MIN_POCET_HER_WIN_BOB:MAX_POCET_HER_WIN_BOB), k , k-1), k - min(MIN_POCET_HER_WIN_BOB:MAX_POCET_HER_WIN_BOB))
    }                                                                 # vzorec pro scitani pravdepodobnosti vyhry boba pri zadanem stavu BODY_ALICE : BODY_BOB
  }
}

E_Bwin_FINAL = E_Bwin/r                                               # prumerna pravdepodobnost vyhry Boba pri stavu 5:3
E_Awin_FINAL = 1 - E_Bwin_FINAL                                       # prumerna pravdepodobnost vyhry Alice pri stavu 5:3

cat("Prumerna pravdepodobnost vyhry Boba pri stavu", BODY_ALICE ,":", BODY_BOB, "=", E_Bwin_FINAL)
cat("Prumerna pravdepodobnost vyhry Alice pri stavu", BODY_ALICE, ":", BODY_BOB, "=", E_Awin_FINAL)
cat("Ferovy podil sanci A:B pro rozdeleni vyhry za stavu", BODY_ALICE, ":", BODY_BOB, "=", E_Awin_FINAL/E_Bwin_FINAL)
cat("Pocet stavu", BODY_ALICE ,":", BODY_BOB, "=", r)
cat("Relativni zastoupeni poctu stavu", BODY_ALICE,":", BODY_BOB, "=", r/s)





