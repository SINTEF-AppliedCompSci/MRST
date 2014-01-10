<?xml version="1.0"?>
<xsl:stylesheet 
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
	xmlns:svg="urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"
	xmlns:draw="urn:oasis:names:tc:opendocument:xmlns:drawing:1.0"
	version="1.0">
	<xsl:output method="text"/>
	<xsl:template match="draw:line">
		      <xsl:value-of select="substring-after(@draw:style-name,'gr')"/>
		      <xsl:text>&#xa;</xsl:text>
        </xsl:template>
</xsl:stylesheet>