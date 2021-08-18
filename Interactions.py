import pandas as pd 
import numpy as np

def All_Contact(df):
    allcontact = df[df['interaction type'].str.contains('contact')]
    allcontact['total'] = allcontact.sum(axis=1)
    allcontact['normalized'] = (allcontact['total'] /9) 
    return(allcontact)

def Charged(df):
    polarcharged = df[df['interaction type'].str.contains('charged')]
    polarcharged['total'] = polarcharged.sum(axis=1)
    polarcharged['normalized'] = (polarcharged['total'] /9) 
    return(polarcharged)
    
def Polar(df):
    polar = df[df['interaction type'].str.contains('polar')]
    polar['total'] = polar.sum(axis=1)
    polar['normalized'] = (polar['total'] /9) 
    return(polar)
    
def Hydrophobic(df):
    hydrophobic = df[df['interaction type'].str.contains('hydrophobic')]
    hydrophobic['total'] = hydrophobic.sum(axis=1)
    hydrophobic['normalized'] = (hydrophobic['total'] /9)
    return(hydrophobic)
    
def Backbone(df):
    Backbone = df[df['interaction type'].str.contains('backbone')]
    Backbone['total'] = Backbone.sum(axis=1)
    Backbone['normalized'] = (Backbone['total'] /9)
    return(Backbone)
