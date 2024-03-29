import re


class GffRow(object):
    def __init__(self, row):
        (
            self.seqid,
            self.source,
            self.type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            self.attributes,
        ) = row.strip().split("\t")
        self.is_gene = self.type == "gene"
        self.parsed_attributes = self.parse_attributes()

    def parse_attributes(self):
        parsed_attributes = {}
        for attr in self.attributes.strip(";").split(";"):
            mm = re.search(r"(?P<key>.*?)(\s+|=)(?P<value>.*)", attr.strip())
            if mm:
                parsed_attributes[mm.group("key").strip()] = mm.group("value").strip().strip('"')
        return parsed_attributes

    def __getitem__(self, key):
        return self.parsed_attributes[key]
