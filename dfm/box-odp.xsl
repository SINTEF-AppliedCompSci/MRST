<?xml version="1.0"?>
<xsl:stylesheet
   xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
   xmlns:style="urn:oasis:names:tc:opendocument:xmlns:style:1.0"
   xmlns:fo="urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0"
   version="1.0">
  <xsl:output method="text"/>
  <!-- find the name of the master page layout (default is PM1) -->
  <xsl:variable name="master">
    <xsl:value-of select="//style:master-page[@style:name='Default']/@style:page-layout-name"/>
  </xsl:variable>
  <!-- get parameters from the master page layout -->
  <xsl:template match="style:page-layout">
    <xsl:if test="@style:name = $master">
      <xsl:value-of select="substring-before(style:page-layout-properties/@fo:page-width, 'cm')"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="substring-before(style:page-layout-properties/@fo:page-height, 'cm')"/>
      <xsl:text> </xsl:text>
    </xsl:if>
  </xsl:template>
  <!-- drop text elements from output -->
  <xsl:template match="text()"/>
</xsl:stylesheet>
