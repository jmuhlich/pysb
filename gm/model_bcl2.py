import pysb


Monomer('caspase8', ['SUB', 'MCL1', 'active'], {'active': ['no', 'yes']})

Monomer('BIDC', ['CSP8', 'PRO', 'ANTI', 'active'], {'active': ['no', 'yes']})

Monomer('BIDM', ['CSP8', 'PRO', 'ANTI', 'active'], {'active': ['no', 'yes']})

Monomer('BAXC', ['ACTIVATOR', 'ANTI', 'active'], {'active': ['no', 'yes']})

(defmonomer BAXC "Bax documentation"
  CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer (BAXM membrane) "Bax documentation"
  CDIM CTETRA CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer BIM "Bim documentation"
  CPRO CANTI)
(defmonomer (BAK membrane) "Bak documentation"
  CDIM CTETRA CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer BCLXLC "BCL-XL documentation"
  CBH3 CPRO)
(defmonomer (BCLXLM membrane) "BCL-XL documentation"
  CBH3 CPRO)
(defmonomer (BCL2 membrane) "BCL2 documentation"
  CBH3 CPRO)
;; Added state to MCL1 to account for degradation. Default of MCL1 is active (unless specified).
(defmonomer MCL1C "MCL1 documentation" 
  CCSP8 CBH3 CPRO (asite :states (member :active :degraded) :default :active))
(defmonomer (MCL1M membrane) "MCL1 documentation" 
  CCSP8 CBH3 CPRO (asite :states (member :active :degraded) :default :active))
(defmonomer BAD "BAD documentation"
  CANTI)
(defmonomer NOXA "NOXA documentation"
  CANTI)
(defmonomer MULE "MULE documentation"
  CANTI)
(defmonomer caspase3 "caspase 3"
  CCSP8 CMCL1 (asite :states (member :active :inactive) :default :inactive))
(defmonomer CYTC "cytochrome C documentation"
  CPORE)
(defmonomer SMAC "SMAC/DIABLO documentation"
  CPORE)
(defmonomer (BAXMPORE membrane) "hokey pore"
  CMITOP)
(defmonomer (BAKPORE  membrane) "hokey pore"
  CMITOP)
