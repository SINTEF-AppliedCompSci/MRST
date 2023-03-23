<?xml version="1.0"?>
<xsl:stylesheet 
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
	xmlns:svg="urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"
	xmlns:draw="urn:oasis:names:tc:opendocument:xmlns:drawing:1.0"
	version="1.0">
	<xsl:output method="text"/>
	<xsl:template match="draw:line">
		      <xsl:value-of select="substring-before(@svg:x1,'cm')"/>
		      <xsl:text> </xsl:text>
		      <xsl:value-of select="substring-before(@svg:y1,'cm')"/>
		      <xsl:text> </xsl:text>
		      <xsl:value-of select="substring-before(@svg:x2,'cm')"/>
		      <xsl:text> </xsl:text>
		      <xsl:value-of select="substring-before(@svg:y2,'cm')"/>
		      <xsl:text>&#xa;</xsl:text>
        </xsl:template>
</xsl:stylesheet>
