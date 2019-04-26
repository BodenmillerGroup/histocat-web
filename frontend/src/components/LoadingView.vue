<template>
  <v-container fluid fill-height>
      <v-layout column align-center justify-center>
        <v-progress-circular
          indeterminate
          color="primary"
        />
        <h5>{{text}}</h5>
      </v-layout>
    </v-container>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';

  @Component
  export default class LoadingView extends Vue {
    @Prop(String) text: string | undefined;
  }
</script>
