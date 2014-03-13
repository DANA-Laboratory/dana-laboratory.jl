---
layout: post
category : code_comment
module: idealGas.jl
tagline: "thermodynamics calculation with julia"
tags : [julia, thermodynamics, sample]
filename:  CpIdeal
---
{% include JB/setup %}

{% capture source_url %}{{page.module}}/{{page.filename}}.jl{% endcapture %}
{% capture comment_url %}{{page.module}}/{{page.filename}}.txt{% endcapture %}

<p align="right">
{% include {{comment_url}} %}
</p> 
{% highlight julia %}
{% include {{source_url}} %}
{% endhighlight %}