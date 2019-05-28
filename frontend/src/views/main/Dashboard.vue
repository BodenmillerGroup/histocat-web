<template>
  <masonry
    :cols="3"
    :gutter="30"
  >
    <ExperimentCard
      v-for="experiment in experiments"
      :experiment="experiment"
      :key="experiment.id"
    />
  </masonry>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readAdminExperiments } from '@/modules/experiment/getters';
  import { dispatchGetExperiments } from '@/modules/experiment/actions';
  import ExperimentCard from '@/components/ExperimentCard.vue';

  @Component({
    components: { ExperimentCard },
  })
  export default class Dashboard extends Vue {

    get experiments() {
      return readAdminExperiments(this.$store);
    }

    async mounted() {
      await dispatchGetExperiments(this.$store);
    }
  }
</script>
