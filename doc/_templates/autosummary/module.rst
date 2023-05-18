{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. topic:: {{ _('Module Attributes') }}

      .. autosummary::
         :nosignatures:
      {% for item in attributes %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. topic:: {{ _('Functions') }}

      .. autosummary::
         :nosignatures:
      {% for item in functions %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. topic:: {{ _('Classes') }}

      .. autosummary::
         :nosignatures:
      {% for item in classes %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. topic:: {{ _('Exceptions') }}

      .. autosummary::
         :nosignatures:
      {% for item in exceptions %}
         {{ item }}
      {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. topic:: Modules

   .. autosummary::
      :nosignatures:
      :toctree:
      :recursive:
   {% for item in modules %}
      {{ item }}
   {%- endfor %}
{% endif %}
{% endblock %}
