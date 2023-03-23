<?xml version="1.0"?>
<xsl:stylesheet 
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:svg="urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"
	xmlns:style="urn:oasis:names:tc:opendocument:xmlns:style:1.0"
    version="1.0">
	<xsl:output method="text"/>
	<xsl:template match="style:graphic-properties">
		      <xsl:value-of select="substring-after(@svg:stroke-color,'#')"/>
		      <xsl:text> </xsl:text>
        </xsl:template>
</xsl:stylesheet>
