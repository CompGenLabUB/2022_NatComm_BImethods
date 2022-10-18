#!/usr/bin/python3
from intermine.webservice import Service
service = Service("https://planmine.mpibpc.mpg.de/planmine/service")
query = service.new_query("Transcript")
query.add_view(
    "primaryIdentifier",
    "symbol",
    "ontologyAnnotations.ontologyTerm.description",
    "ontologyAnnotations.ontologyTerm.identifier",
    "ontologyAnnotations.ontologyTerm.name",
    "ontologyAnnotations.ontologyTerm.namespace",
    "ontologyAnnotations.ontologyTerm.obsolete"
)
query.add_sort_order("primaryIdentifier", "ASC")

for row in query.rows():
    print(row["primaryIdentifier"], \
          row["symbol"], \
          row["ontologyAnnotations.ontologyTerm.namespace"], \
          row["ontologyAnnotations.ontologyTerm.identifier"], \
          row["ontologyAnnotations.ontologyTerm.name"], \
          row["ontologyAnnotations.ontologyTerm.description"], \
          row["ontologyAnnotations.ontologyTerm.obsolete"],
          sep='\t')
