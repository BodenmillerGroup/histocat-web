<template>
  <v-card class="ma-3 pa-3" tile hover>
    <v-card-title>
      <v-layout column>
        <h5 class="headline">{{experiment.name}}</h5>
        <span class="caption"><v-icon small>mdi-calendar-outline</v-icon> {{createdAt}}</span>
      </v-layout>
    </v-card-title>
    <v-card-text>
      {{experiment.description}}
    </v-card-text>
    <v-card-actions>
      <v-btn :to="{name: 'main-experiment', params: {id: experiment.id}}">View</v-btn>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';
  import { IExperiment } from '@/modules/experiment/models';

  @Component
  export default class ExperimentCard extends Vue {
    @Prop(Object) experiment?: IExperiment;

    get createdAt() {
      return this.experiment && new Date(this.experiment.created_at).toLocaleString('de-ch', { timeZone: 'UTC' });
    }
  }
</script>
