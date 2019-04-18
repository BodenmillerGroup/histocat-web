<template>
  <v-container fluid>
    <v-card class="ma-3 pa-3">
      <v-card-title primary-title>
        <div class="headline primary--text">Edit Experiment</div>
      </v-card-title>
      <v-card-text>
        <template>
          <v-form
            v-model="valid"
            ref="form"
            lazy-validation
          >
            <v-text-field
              label="Name"
              v-model="name"
              v-validate="'required'"
              data-vv-name="name"
              :error-messages="errors.collect('name')"
              required
            ></v-text-field>
            <v-text-field
              label="Description"
              v-model="description"
            ></v-text-field>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn
          @click="submit"
          :disabled="!valid"
        >
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { IExperimentUpdate } from '@/modules/experiment/models';
  import { dispatchGetExperiments, dispatchUpdateExperiment } from '@/modules/experiment/actions';
  import { readAdminOneExperiment } from '@/modules/experiment/getters';

  @Component
  export default class EditExperiment extends Vue {
    valid = true;
    name: string = '';
    description: string = '';

    async mounted() {
      await dispatchGetExperiments(this.$store);
      this.reset();
    }

    reset() {
      this.name = '';
      this.description = '';
      this.$validator.reset();
      if (this.experiment) {
        this.name = this.experiment.name;
        this.description = this.experiment.description;
      }
    }

    cancel() {
      this.$router.back();
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        const data: IExperimentUpdate = {};
        if (this.name) {
          data.name = this.name;
        }
        if (this.description) {
          data.description = this.description;
        }
        await dispatchUpdateExperiment(this.$store, { id: this.experiment!.id, data: data });
        this.$router.push('/main/admin/experiments');
      }
    }

    get experiment() {
      return readAdminOneExperiment(this.$store)(+this.$router.currentRoute.params.id);
    }
  }
</script>
