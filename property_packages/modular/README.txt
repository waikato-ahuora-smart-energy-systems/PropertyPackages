Module Version 0.1


--------------------------------------------------
Purpose of this module is to generate a json file
parsable by idaes GenericParameterBlock.
--------------------------------------------------

Design Decisions

--------------------------------------------------
1 - plug and play design
--------------------------------------------------
The design should be created in a way that allows
future configuration files to be implemented despite
the drastic differences in structure
--------------------------------------------------
- Idea
--------------------------------------------------

--------------------------------------------------
Layout Files (Multiple Property Packages)
--------------------------------------------------

Layout Files
| - Peng-Robinson
| - NRTL

Within the layout files we have a system similar to the
GenericParameterBlock pattern. For example using peng-robinson

"PR": {
  "base_units": BuildMethod
  "parameter_data": BuildMethod
}

Structure of Build-Methods

Accepts -> list of compounds
Returns -> json file structure

Layout Files Have a Link to their corresponding property package

Build methods can be organised in the following ways

Build Methods
- base
- | - SI_base_units
- | - REF_TEMP
- salts
- electrolyte

Custom build methods?

- What happens when we need to generate -> parameter_data
- Custom Logic?

--------------------------------------------------
Retrieving information
--------------------------------------------------
Information needs to be retrieved for individual
compounds. The individual compound information can
be found in the data files.

1- Compound names passed in a validated
2- Compound names are converted into a list of Compound objects

alternative: compounds are classes

--------------------------------------------------
Compound Verification
--------------------------------------------------



Key Points to ask about
--------------------------------------------------
- best method of storing compound information for retrieval and use.
- best way of actually tracking what compounds are available (for input validation)




