"""
Function:  queryActivities 

    gets information for all relevant actives from ChEMBL
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""  

def activities(species):

    query ="SELECT act.molregno                                               \
            ,act.standard_value                                               \
            ,act.standard_type                                                \
            ,act.standard_units                                               \
            ,ass.assay_type                                                   \
            ,ass.assay_id                                                     \
            ,td.pref_name                                                     \
            ,td.protein_accession                                             \
            ,docs.doc_id                                                      \
                                                                              \
      FROM assay2target a2t                                                   \
        JOIN (SELECT a2t.assay_id                                             \
                FROM assay2target a2t                                         \
               GROUP BY assay_id HAVING COUNT(a2t.tid) =1) selass             \
            ON selass.assay_id = a2t.assay_id                                 \
         JOIN activities act                                                  \
            ON a2t.assay_id = act.assay_id                                    \
          JOIN target_dictionary td                                           \
            ON a2t.tid = td.tid                                               \
          JOIN assays ass                                                     \
            ON ass.assay_id = act.assay_id                                    \
          JOIN docs                                                           \
            ON act.doc_id = docs.doc_id                                       \
                                                                              \
        WHERE td.target_type = 'PROTEIN'                                      \
          AND ass.assay_type='B'                                              \
          AND act.relation ='='                                               \
          AND a2t.multi=0                                                     \
          AND a2t.complex=0                                                   \
          AND a2t.relationship_type = 'D'                                     \
          AND act.standard_type IN('Ki','IC50', 'EC50', 'pA2','pKi')          \
          AND organism = '%s'" % species
    
    return query
 #

def paralogs(paralogTable):
    inputFile = open(paralogTable, 'r') 
    paralogDict={}
    lines = inputFile.readlines()
    for line in lines:
        elements = line.split('\t')
        paralogDict[elements[0]]=0
        paralogDict[elements[1]]=0
    paralogString = "\',\'".join(paralogDict.keys())
    query = "SELECT act.molregno                                              \
            ,act.standard_value                                               \
            ,act.standard_type                                                \
            ,act.standard_units                                               \
            ,ass.assay_type                                                   \
            ,ass.assay_id                                                     \
            ,td.pref_name                                                     \
            ,td.protein_accession                                             \
            ,docs.doc_id                                                      \
        FROM assay2target a2t                                                 \
          JOIN (SELECT a2t.assay_id                                           \
                FROM assay2target a2t                                         \
               GROUP BY assay_id HAVING COUNT(a2t.tid) =1) selass             \
            ON selass.assay_id = a2t.assay_id                                 \
          JOIN activities act                                                 \
            ON a2t.assay_id = act.assay_id                                    \
          JOIN target_dictionary td                                           \
            ON a2t.tid = td.tid                                               \
          JOIN assays ass                                                     \
            ON ass.assay_id = act.assay_id                                    \
          JOIN docs                                                           \
            ON act.doc_id = docs.doc_id                                       \
        WHERE td.target_type = 'PROTEIN'                                      \
          AND ass.assay_type='B'                                              \
          AND act.relation ='='                                               \
          AND a2t.multi=0                                                     \
          AND a2t.complex=0                                                   \
          AND a2t.relationship_type = 'D'                                     \
          AND act.standard_type IN('Ki','IC50', 'EC50', 'pA2','pKi')          \
          AND td.protein_accession IN('%s')"%paralogString
    
    return query
	

