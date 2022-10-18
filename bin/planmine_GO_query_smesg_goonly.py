#!/usr/bin/python3
from intermine.webservice import Service
service = Service("https://planmine.mpibpc.mpg.de/planmine/service")
query = service.new_query("PredictedGene")
query.add_view(
    "primaryIdentifier",
    "symbol",
    "transcripts.ontologyAnnotations.ontologyTerm.description",
    "transcripts.ontologyAnnotations.ontologyTerm.identifier",
    "transcripts.ontologyAnnotations.ontologyTerm.name",
    "transcripts.ontologyAnnotations.ontologyTerm.namespace",
    "transcripts.ontologyAnnotations.ontologyTerm.obsolete"
)
query.add_sort_order("primaryIdentifier", "ASC")

for row in query.rows():
    print(row["primaryIdentifier"], \
          row["symbol"], \
          row["transcripts.ontologyAnnotations.ontologyTerm.namespace"], \
          row["transcripts.ontologyAnnotations.ontologyTerm.identifier"], \
          row["transcripts.ontologyAnnotations.ontologyTerm.name"], \
          row["transcripts.ontologyAnnotations.ontologyTerm.description"], \
          row["transcripts.ontologyAnnotations.ontologyTerm.obsolete"],
          sep='\t')
