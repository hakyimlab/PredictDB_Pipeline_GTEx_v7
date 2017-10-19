import sys
from sqlalchemy import create_engine, Column, Decimal, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()

class ModelSummary(Base):
    __tablename__ = 'model_summaries'

    gene = Column(String, primary_key=True)
    gene_name = Column(String(20))
    gene_type = Column(String(10))
    chromosome = Column(String(2))
    start_pos = Column(Integer)
    end_pos = Column(Integer)

    alpha = Column(Integer)
    
class SNPInfo(Base):
    __tablename__ = 'snp_info'

    rsid = Column(String, primary_key=True)
    chromosome = Column(String(2))
    position = Column(Integer)
    reference_allele = Column(String(1))
    effect_allele = Column(String(1))
    minor_allele_frequency = Column(Decimal)
    
    rsid_cov1 = relationship('Covariance', 

class Weight(Base):
    __tablename__ = 'weights'

    gene = Column(String, ForeignKey('model_summaries.gene'))
    rsid = Column(String, ForeignKey('snp_info.rsid'))
    beta = Column(Decimal)

class Covariance(Base):
    __tablename__ = 'covariances'
    
    gene = Column(String, ForeignKey('model_summaries.gene'), primary_key=True)
    rsid1 = Column(String, ForeignKey('snp_info.rsid'), primary_key=True)
    rsid2 = Column(String, ForeignKey('snp_info.rsid'), primary_key=True)
    covariance = Column(Decimal)

    model_summary = relationship("ModelSummary")
    

class ConstructionInfo(Base):
    __tablename__ = 'construction'

    chromosome = Column(String(2), primary_key=True)
    seed = Column(Integer)

class SampleInfo(Base):
    __tablename__ = 'sample_info'
    
    n_samples = Column(Integer)
    population = Column(String(30))
    tissue_sample = Column(String(50), primary_key=True)

